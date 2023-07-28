//compile with OpenCilk-10.0.1-Linux/bin/clang -fopencilk -pthread task_1_opencilk.c mmio.c -o task_1_opencilk
//and CILK_NWORKERS=XX time ./task_1_opencilk

#include "cilk/cilk.h"
#include "mmio.h"
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define graphName "mycielskian13.mtx"
#define BILLION 1E9
#define CLOCK_REALTIME 0

pthread_mutex_t myMutex;

/*Triangle : function that takes as input a matrix A in CSC format 
and give as output a matrix C that is the product of A.*(A*A) also
 in CSC format .Additionally it computes the number val of triangles 
 between the nodes of the graph that matrix A represents (val = sum(elements of c) 
 where c=C*e (e :vector of length n with ones)) We now use Opencilk add-on to
 help us with the parallelization  */

void triangle(int* rowind, int* colptr, int N, int nz, int* cCol, int* cRow, int* cVal, int* valPtr)
{
    int val = 0;
    int currentVal;
    cCol[0] = 0;

    int currentPtr;
    int nextPtr;
    int col;
    int row;

    for (int i = 0; i < N; i++)

    {
        currentPtr = colptr[i];
        nextPtr = colptr[i + 1];

        for (int j = currentPtr; j < nextPtr; j++) {
            col = i;
            row = rowind[j];

            cRow[colptr[i] + (j - currentPtr)] = row;

            currentVal = 0;
            cilk_for(int h = colptr[col]; h < colptr[col + 1]; h++)
            {
                for (int k = colptr[row]; k < colptr[row + 1]; k++) {

                    if (rowind[h] == rowind[k]) {
                        pthread_mutex_lock(&myMutex);
                        currentVal++;
                        (*valPtr)++;
                        pthread_mutex_unlock(&myMutex);
                    }
                }
            }
            cVal[colptr[i] + (j - currentPtr)] = currentVal;
        }
        cCol[i + 1] = colptr[i + 1];
    }

    printf("C matrix in CSC format created\n");
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

    //Start the clock

    struct timespec requestStart, requestEnd;
    clock_gettime(CLOCK_REALTIME, &requestStart);

    //Function triangle

    triangle(csc_row, csc_col, N, nz, cCol, cRow, cVal, valPtr);
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
    /* for (int i = 0; i < nz; i++)

    {
        if (i < N) {
            printf("c %d r %d v %d\n", cCol[i], cRow[i], cVal[i]);
        } else {
            printf(" r %d v %d\n", cRow[i], cVal[i]);
        }
    }*/

    return 0;
}
