//compile with gcc -Wall -g p_thr.c mmio.c -o p_thr -pthread

#include "mmio.h"
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define graphName "dblp-2010.mtx"
#define BILLION 1E9
#define NUMOFTHREADS 75
pthread_mutex_t myMutex;
//Define a structure for the triangle function

typedef struct argumentsOfTriangle {
    int* rowind;
    int* colptr;
    int N;
    int nz;
    int* cCol;
    int* cRow;
    int* cVal;
    int* valPtr;
    int* Npartitions;
    int threadID;
};

/*Triangle : function that takes as input a matrix A in CSC format 
and give as output a matrix C that is the product of A.*(A*A) also
 in CSC format .Additionally it computes the number val of triangles 
 between the nodes of the graph that matrix A represents (val = sum(elements of c) 
 where c=C*e (e :vector of length n with ones)) Here we have adapt the function to be suitable for pthreads
 So we pass the arguments through a struct Because of the kind of the algorith we only have to 
 split the first loop in as many parts as the threads  */

void* triangle(void* received_Struct)
{
    struct argumentsOfTriangle* argofT = (struct argumentsOfTriangle*)received_Struct;

    int currentVal;
    argofT->cCol[0] = 0;

    int currentPtr;
    int nextPtr;
    int col;
    int row;

    for (int i = argofT->Npartitions[argofT->threadID]; i < argofT->Npartitions[argofT->threadID + 1]; i++) {
        currentPtr = argofT->colptr[i];
        nextPtr = argofT->colptr[i + 1];

        for (int j = currentPtr; j < nextPtr; j++) {
            col = i;
            row = argofT->rowind[j];

            argofT->cRow[argofT->colptr[i] + (j - currentPtr)] = row;

            currentVal = 0;
            for (int h = argofT->colptr[col]; h < argofT->colptr[col + 1]; h++) {
                for (int k = argofT->colptr[row]; k < argofT->colptr[row + 1]; k++) {

                    if (argofT->rowind[h] == argofT->rowind[k]) {

                        currentVal++;
                        pthread_mutex_lock(&myMutex);   //we use mutex to confront race events
                        (*argofT->valPtr)++;
                        pthread_mutex_unlock(&myMutex);
                    }
                }
            }
            argofT->cVal[argofT->colptr[i] + (j - currentPtr)] = currentVal;
        }
        argofT->cCol[i + 1] = argofT->colptr[i + 1];
    }

    // printf("C matrix in CSC format created , Thread: %d \n", argofT->threadID);
}

//Function from github that transform a COO matrix to a CSC

void coo2csc(
    uint32_t* const row, /*!< CSC row start indices */
    uint32_t* const col, /*!< CSC column indices */
    uint32_t const* const row_coo, /*!< COO row indices */
    uint32_t const* const col_coo, /*!< COO column indices */
    uint32_t const nnz, /*!< Number of nonzero elements */
    uint32_t const n, /*!< Number of rows/columns */
    uint32_t const isOneBased /*!< Whether COO is 0- or 1-based */
)
{

    // ----- cannot assume that input is already 0!
    for (uint32_t l = 0; l < n + 1; l++)
        col[l] = 0;

    // ----- find the correct column sizes
    for (uint32_t l = 0; l < nnz; l++)
        col[col_coo[l] - isOneBased]++;

    // ----- cumulative sum
    for (uint32_t i = 0, cumsum = 0; i < n; i++) {
        uint32_t temp = col[i];
        col[i] = cumsum;
        cumsum += temp;
    }
    col[n] = nnz;
    // ----- copy the row indices to the correct place
    for (uint32_t l = 0; l < nnz; l++) {
        uint32_t col_l;
        col_l = col_coo[l] - isOneBased;

        uint32_t dst = col[col_l];
        row[dst] = row_coo[l] - isOneBased;

        col[col_l]++;
    }
    // ----- revert the column pointers
    for (uint32_t i = 0, last = 0; i < n; i++) {
        uint32_t temp = col[i];
        col[i] = last;
        last = temp;
    }
}

int main(int argc, char* argv[])
{
    int val = 0;
    int* valPtr;
    int ret_code;
    int M, N, nz;
    int i, *I, *J;
    int* Npartitions;
    FILE* f;
    MM_typecode matcode;

    int *cCol, *cRow, *cVal;

    valPtr = &val;
    //reading the graph and save the informations

    if (fopen(graphName, "r") == NULL) {
        printf("Failure in reading file\n");
    }
    f = fopen(graphName, "r");

    int bannerOutput;
    bannerOutput = mm_read_banner(f, &matcode);

    if (bannerOutput == 1) {
        printf("Successfull Banner \n");
    }
    if (mm_is_pattern(matcode) == 0) {
        printf("Pattern? %d \n", mm_is_pattern(matcode)); //Check if the matrix is Pattern
    }
    ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz);

    I = (int*)malloc(nz * sizeof(int));
    J = (int*)malloc(nz * sizeof(int));

    for (i = 0; i < nz; i++) {
        fscanf(f, "%d %d \n", &I[i], &J[i]);
        I[i]--; // adjust from 1-based to 0-based
        J[i]--;
    }

    if (f != stdin)
        fclose(f);
    // coo to csc

    uint32_t* csc_row = (uint32_t*)malloc(nz * sizeof(uint32_t));
    uint32_t* csc_col = (uint32_t*)malloc((N + 1) * sizeof(uint32_t));

    coo2csc(csc_row, csc_col, J, I, nz, N, 0);

    //Creating the matrices for C CSC format save

    cCol = (int*)malloc((N + 1) * sizeof(int));
    cRow = (int*)malloc(nz * sizeof(int));
    cVal = (int*)malloc(nz * sizeof(int));

    //Create the repeat numbers for every thread

    int range = floor(N / NUMOFTHREADS);
    Npartitions = (int*)malloc((NUMOFTHREADS + 1) * sizeof(int));
    Npartitions[0] = 0;
    for (int i = 1; i < NUMOFTHREADS; i++) {
        Npartitions[i] = Npartitions[i - 1] + range;
    }
    Npartitions[NUMOFTHREADS] = N;
    /*
    for (int i = 0; i < NUMOFTHREADS + 1; i++) {
        printf("%d\n", Npartitions[i]);
    }
*/
    //Seting up the threads
    int* ids;
    ids = (int*)malloc(NUMOFTHREADS * sizeof(int));

    pthread_t* threads;
    threads = (pthread_t*)malloc(NUMOFTHREADS * sizeof(pthread_t));

    struct argumentsOfTriangle* matrixOfargumentsOfTriangle;
    matrixOfargumentsOfTriangle = malloc(NUMOFTHREADS * sizeof(struct argumentsOfTriangle));

    pthread_attr_t pthread_custom_attr;
    pthread_attr_init(&pthread_custom_attr);

    //Define a structure for the triangle function

    struct argumentsOfTriangle argT1;
    argT1.rowind = csc_row;
    argT1.colptr = csc_col;
    argT1.N = N;
    argT1.nz = nz;
    argT1.cCol = cCol;
    argT1.cRow = cRow;
    argT1.cVal = cVal;
    argT1.valPtr = valPtr;
    argT1.Npartitions = Npartitions;

    //Start the clock

    struct timespec requestStart, requestEnd;
    clock_gettime(CLOCK_REALTIME, &requestStart);

    //Function triangle

    for (int i = 0; i < NUMOFTHREADS; i++) {
        ids[i] = i;
        matrixOfargumentsOfTriangle[i] = argT1;
        matrixOfargumentsOfTriangle[i].threadID = ids[i];
        pthread_create(&threads[i], &pthread_custom_attr, triangle, &matrixOfargumentsOfTriangle[i]);
    }
    for (i = 0; i < NUMOFTHREADS; i++) {
        pthread_join(threads[i], NULL);
    }

    printf("Number of Triangles %d\n", val);
    /*The real number of triangles is 6 times smallerbecause we 
    count the same triangle with different orientation 6 different times*/

    // Calculate time it took

    clock_gettime(CLOCK_REALTIME, &requestEnd);

    double accum = (requestEnd.tv_sec - requestStart.tv_sec)
        + (requestEnd.tv_nsec - requestStart.tv_nsec)
            / BILLION;
    printf("Time it took :%lf seconds\n", accum);

    //The code in comments prints the 3 matrices of C CSC format
    /*for (int i = 0; i < nz; i++)

    {
        if (i < N) {
            printf("c %d r %d v %d\n", cCol[i], cRow[i], cVal[i]);
        } else {
            printf(" r %d v %d\n", cRow[i], cVal[i]);
        }
    }*/

    return 0;
}
