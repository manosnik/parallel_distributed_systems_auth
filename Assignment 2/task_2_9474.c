//compile with mpicc task_2_9474.c -o task_2_9474
//compile with mpirun -np x ./task_2_9474
#include "mpich/mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

//How many datas and with how many columns (=vector's dimensions)
#define ROWS 8192
#define COL 15
#define FILENAME "Dry_Bean_Dataset2.txt"
#define BILLION 1E9

/*Function that computes the euclidian distance ^2 between two vectors
of d dimensions */

double euclidianDistance(double* center, double* point, int d)
{
    double result = 0;
    for (int i = 0; i < d; i++) {
        result = result + (point[i] - center[i]) * (point[i] - center[i]);
    }
    return result;
}

/*The function swap ,partition ,quickselect ,ksmallest are all contribute 
in finding the k th smallest element in an array through the quickselect 
sort algorithm*/

void swap(double* a, double* b)
{
    double temp = *a;
    *a = *b;
    *b = temp;
}

int partition(double* A, int left, int right)
{
    double pivot = A[right];
    int i = left, x;
    for (x = left; x < right; x++) {
        if (A[x] < pivot) {
            swap(&A[i], &A[x]);
            i++;
        }
    }

    swap(&A[i], &A[right]);
    return i;
}

double quickselect(double* A, int left, int right, int k)
{

    int p = partition(A, left, right);

    //k equals pivot got lucky
    if (p == k - 1) {
        return A[p];
    }
    //k less than pivot
    else if (k - 1 < p) {
        return quickselect(A, left, p - 1, k);
    }
    //k greater than pivot
    else {
        return quickselect(A, p + 1, right, k);
    }
}

double ksmallest(double* A, int n, int k)
{
    int left = 0;
    int right = n - 1;

    return quickselect(A, left, right, k);
}

/*It computes 2 arrays
1. The 3D finalArr that the finalArr[i][j][k] (k==0 or k==1) element means 
that the i th rank has to send (and receive)  finalArr[i][j][1] elements
to finalArr[i][j][0] process.
2. countArr that consists of how many exchanges has the i th process to execute*/

void function(int* exchangeTable, int world_size, int* countArr, int*** finalArr)
{

    for (int i = 0; i < world_size; i++) {
        countArr[i] = 0;
    }
    int* exchangeTableCopy = (int*)malloc(world_size * sizeof(int));

    for (int i = 0; i < world_size; i++) {
        exchangeTableCopy[i] = exchangeTable[i];
    }

    int count = 0;
    for (int i = 0; i < world_size / 2; i++) {
        count = 0;

        for (int j = world_size / 2; j < world_size; j++) {
            if (exchangeTable[i] > 0) {
                if (exchangeTable[i] <= exchangeTable[j]) {
                    exchangeTable[j] = exchangeTable[j] - exchangeTable[i];
                    exchangeTable[i] = 0;
                    countArr[j]++;
                    countArr[i]++;
                } else {
                    if (exchangeTable[j] > 0) {
                        exchangeTable[i] = exchangeTable[i] - exchangeTable[j];
                        exchangeTable[j] = 0;
                        countArr[i]++;
                        countArr[j]++;
                    }
                }
            }
        }
    }

    int* countersOfSecondHalf = (int*)malloc(world_size * sizeof(int));
    for (int i = 0; i < world_size; i++) {
        countersOfSecondHalf[i] = 0;
    }

    for (int i = 0; i < world_size; i++) {
        finalArr[i] = (int**)malloc(countArr[i] * sizeof(int*));
        for (int k = 0; k < countArr[i]; k++) {
            finalArr[i][k] = (int*)malloc(2 * sizeof(int));
        }
    }
    for (int i = 0; i < world_size; i++) {

        if (i < world_size / 2) {

            count = 0;
            for (int j = world_size / 2; j < world_size; j++) {
                if (exchangeTableCopy[i] > 0) {
                    if (exchangeTableCopy[i] <= exchangeTableCopy[j]) {

                        finalArr[i][count][0] = j;
                        finalArr[i][count][1] = exchangeTableCopy[i];
                        finalArr[j][countersOfSecondHalf[j]][0] = i;
                        finalArr[j][countersOfSecondHalf[j]][1] = exchangeTableCopy[i];
                        exchangeTableCopy[j] = exchangeTableCopy[j] - exchangeTableCopy[i];
                        exchangeTableCopy[i] = 0;
                        countersOfSecondHalf[j]++;
                        count++;
                    } else {
                        if (exchangeTableCopy[j] > 0) {
                            finalArr[i][count][0] = j;
                            finalArr[i][count][1] = exchangeTableCopy[j];
                            finalArr[j][countersOfSecondHalf[j]][0] = i;
                            finalArr[j][countersOfSecondHalf[j]][1] = exchangeTableCopy[j];
                            exchangeTableCopy[i] = exchangeTableCopy[i] - exchangeTableCopy[j];
                            exchangeTableCopy[j] = 0;
                            countersOfSecondHalf[j]++;
                            count++;
                        }
                    }
                }
            }
        }
    }
}

/*Recursive function that redestributes the datas in order that the p/2 first processes 
has the n/2 elements with the smallest distance from pivot and the p/2 last the  n/2 
elements with the bigest distance from pivot*/

void distributeByMedian(int chiefRank, MPI_Comm mpiCommunicator, int numberOfData, int numberOfCol, int numberOfRows, double dataPortion[numberOfRows][numberOfCol], double* pivot)
{

    MPI_Barrier(mpiCommunicator);

    int world_size;
    int world_rank;
    MPI_Comm_size(mpiCommunicator, &world_size);
    MPI_Comm_rank(mpiCommunicator, &world_rank);

    //if the processes that get involved are less than 2 the recursion is over

    if (world_size > 1) {

        //choosing the pivot - only the first time- Does not matter which element is selected

        double quickselectOutput;
        double* allDistances;
        if (world_rank == chiefRank) {
            if (mpiCommunicator == MPI_COMM_WORLD) {
                //  srand(time(0));
                for (int i = 0; i < numberOfCol; i++) {
                    pivot[i] = dataPortion[3][i]; //rand() % (numberOfData / world_size)
                }
            }
            allDistances = (double*)malloc((numberOfData * sizeof(double)));
        }
        //Sending the pivot to everyone and they compute their distances

        MPI_Bcast(pivot, numberOfCol, MPI_DOUBLE, chiefRank, mpiCommunicator);
        double* distances;
        distances = (double*)malloc((numberOfData / world_size) * sizeof(double));
        for (int i = 0; i < (numberOfData / world_size); i++) {
            distances[i] = euclidianDistance(pivot, dataPortion[i], numberOfCol);
        }
        MPI_Barrier(mpiCommunicator);

        //The chief collects all the distances , gather them in an array , computes the median value and sends it to everyone

        MPI_Gather(distances, (numberOfData / world_size), MPI_DOUBLE, allDistances, (numberOfData / world_size), MPI_DOUBLE, chiefRank, mpiCommunicator);

        if (world_rank == chiefRank) {
            quickselectOutput = ksmallest(allDistances, numberOfData, (int)(numberOfData / 2));

            for (int i = 0; i < numberOfData; i++) {
                if (allDistances[i] == quickselectOutput && i != (numberOfData / 2) - 1) { // problem if there are identical datas
                    printf("Problem with your Database i (%d)\n", i);
                }
            }
        }

        MPI_Barrier(mpiCommunicator);

        MPI_Bcast(&quickselectOutput, 1, MPI_DOUBLE, chiefRank, mpiCommunicator);

        //Every process finds how many of its distances are lower than pivot distance and how many are higher
        //Then with PreFixScan we put the datas with higher distance in the first sector of the array

        int lowerCounter = 0, upperCounter = 0, counter = 0;

        int* y;
        y = (int*)malloc((numberOfData / world_size) * sizeof(int));

        for (int i = 0; i < (numberOfData / world_size); i++) {

            if (quickselectOutput >= distances[i]) {
                lowerCounter++;
                y[i] = 0;
            }
            if (quickselectOutput < distances[i]) {
                upperCounter++;
                y[i] = 1;
            }
        }
        int* z;
        z = (int*)malloc((numberOfData / world_size) * sizeof(int));

        for (int i = 0; i < (numberOfData / world_size); i++) {
            z[i] = 0;
            for (int j = 0; j <= i; j++) {
                if (y[i] == 1 && y[j] == 1) {
                    z[i]++;
                }
            }
        }

        for (int i = 0; i < (numberOfData / world_size); i++) {
            if (y[i] == 1) {
                swap(&distances[i], &distances[z[i] - 1]);
                for (int j = 0; j < numberOfCol; j++) {
                    swap(&dataPortion[i][j], &dataPortion[z[i] - 1][j]);
                }
            }
        }

        //Every process finds out how many elements it has to swap

        int* exchangeTable;
        int exchangeValue;
        if (world_rank == chiefRank) {
            exchangeTable = (int*)malloc(world_size * sizeof(int));
            if (chiefRank < world_size / 2) {
                exchangeValue = upperCounter;
            } else {
                exchangeValue = lowerCounter;
            }
        } else {
            if (world_rank < world_size / 2) {
                exchangeValue = upperCounter;
            } else {
                exchangeValue = lowerCounter;
            }
        }

        //Initial creation of the exchange Table
        MPI_Barrier(mpiCommunicator);
        MPI_Gather(&exchangeValue, 1, MPI_INT, exchangeTable, 1, MPI_INT, chiefRank, mpiCommunicator);

        // With the help of the function function we make a 3D exchange array

        int* countArr;
        int*** finalArr;
        if (world_rank == chiefRank) {
            countArr = (int*)malloc(world_size * sizeof(int));
            finalArr = (int***)malloc(world_size * sizeof(int**));
            function(exchangeTable, world_size, countArr, finalArr);
        }

        MPI_Barrier(mpiCommunicator);

        //The chiefRank informs the other about how many sends/recieves they will have to make

        int personalCountArr;
        MPI_Scatter(countArr, 1, MPI_INT, &personalCountArr, 1, MPI_INT, chiefRank, mpiCommunicator);
        MPI_Barrier(mpiCommunicator);

        // Making the 3D array 1D
        // Computes how many memory positions they will have to exchange and informs every process

        int* finalArr1D;
        int* displ;
        int* countArrX2;
        if (world_rank == chiefRank) {
            int sumOfCountArr = 0;
            for (int i = 0; i < world_size; i++) {
                sumOfCountArr = sumOfCountArr + countArr[i];
            }

            displ = (int*)malloc(world_size * sizeof(int));
            finalArr1D = (int*)malloc(sumOfCountArr * 2 * sizeof(int));
            int r = 0;
            for (int i = 0; i < world_size; i++) {
                displ[i] = r;
                for (int k = 0; k < countArr[i]; k++) {
                    for (int t = 0; t < 2; t++) {
                        finalArr1D[r] = finalArr[i][k][t];
                        r++;
                    }
                }
            }

            countArrX2 = (int*)malloc(world_size * sizeof(int));

            for (int i = 0; i < world_size; i++) {
                countArrX2[i] = countArr[i] * 2;
            }
        }

        int* personalFinalArr1D = (int*)malloc(personalCountArr * 2 * sizeof(int));
        MPI_Scatterv(finalArr1D, countArrX2, displ, MPI_INT, personalFinalArr1D, personalCountArr * 2, MPI_INT, chiefRank, mpiCommunicator);
        MPI_Barrier(mpiCommunicator);

        //Final transfer of the Data
        int numberOfTemporaryDataArrRows;
        if (world_rank < world_size / 2) {
            numberOfTemporaryDataArrRows = upperCounter;

        } else {
            numberOfTemporaryDataArrRows = lowerCounter;
        }

        double temporaryData[numberOfTemporaryDataArrRows][numberOfCol];

        MPI_Request mprequest;
        int sumOfPreviousSentData = 0;
        int sumOfPreviousRecivedData = 0;

        for (int i = 0; i < personalCountArr; i++) {
            if (world_rank < world_size / 2) {
                MPI_Isend(dataPortion[sumOfPreviousSentData], personalFinalArr1D[2 * i + 1] * numberOfCol, MPI_DOUBLE, personalFinalArr1D[2 * i], 30, mpiCommunicator, &mprequest);
                sumOfPreviousSentData = sumOfPreviousSentData + personalFinalArr1D[2 * i + 1];
                MPI_Irecv(temporaryData[sumOfPreviousRecivedData], personalFinalArr1D[2 * i + 1] * numberOfCol, MPI_DOUBLE, personalFinalArr1D[2 * i], 30, mpiCommunicator, &mprequest);
                sumOfPreviousRecivedData = sumOfPreviousRecivedData + personalFinalArr1D[2 * i + 1];

            } else {

                MPI_Isend(dataPortion[upperCounter + sumOfPreviousSentData], personalFinalArr1D[2 * i + 1] * numberOfCol, MPI_DOUBLE, personalFinalArr1D[2 * i], 30, mpiCommunicator, &mprequest);
                sumOfPreviousSentData = sumOfPreviousSentData + personalFinalArr1D[2 * i + 1];
                MPI_Irecv(&temporaryData[sumOfPreviousRecivedData], personalFinalArr1D[2 * i + 1] * numberOfCol, MPI_DOUBLE, personalFinalArr1D[2 * i], 30, mpiCommunicator, &mprequest);
                sumOfPreviousRecivedData = sumOfPreviousRecivedData + personalFinalArr1D[2 * i + 1];
            }
        }

        MPI_Barrier(mpiCommunicator);

        if (world_rank < world_size / 2) {
            for (int i = 0; i < upperCounter; i++) {
                for (int j = 0; j < numberOfCol; j++) {
                    dataPortion[i][j] = temporaryData[i][j];
                }
            }

        } else {

            for (int i = 0; i < lowerCounter; i++) {
                for (int j = 0; j < numberOfCol; j++) {
                    dataPortion[upperCounter + i][j] = temporaryData[i][j];
                }
            }
        }

        MPI_Barrier(mpiCommunicator);

        //Set the scene for the recursion through creating new MPI communicators

        int color = world_rank / (world_size / 2);
        MPI_Comm mpiHalfCommunicator;
        MPI_Comm_split(mpiCommunicator, color, world_rank, &mpiHalfCommunicator);

        distributeByMedian(0, mpiHalfCommunicator, numberOfData / 2, numberOfCol, numberOfRows, dataPortion, pivot);

        MPI_Barrier(mpiCommunicator);
    }
}

int main(int argc, char** argv)
{
    //Start the clock

    struct timespec requestStart, requestEnd;
    clock_gettime(CLOCK_REALTIME, &requestStart);

    //Reading of the data from a  file of ISC database

    int numberOfRows = ROWS;
    int numberOfCol = COL;

    double data[numberOfRows][numberOfCol];

    FILE* dataFile = fopen(FILENAME, "r");
    if (dataFile == NULL) {
        perror("Failure in file reading");
        exit(1);
    }

    char buf[30000];
    char* token;
    char* delimeter = ",";

    for (int i = 0; i < numberOfRows; i++) {
        fgets(buf, 30000, dataFile);
        token = strtok(buf, delimeter);
        data[i][0] = atof(token);
        for (int j = 1; j < numberOfCol; j++) {
            token = strtok(NULL, delimeter);
            data[i][j] = atof(token);
            if (token == NULL)
                break;
        }
    }

    //Some initial informations and the start of MPI

    int numberOfData = numberOfRows;
    int err;
    MPI_Status mpistat;

    double* startDistances = (double*)malloc(numberOfData * sizeof(double));
    int c = 0;

    /*  for (int i = 0; i < numberOfData; i++) {
        if (i % 1000 == 0)
            printf("%d\n", i);
        for (int j = 0; j < numberOfData; j++) {
            startDistances[j] = euclidianDistance(data[i], data[j], COL);
        }

        for (int k = 0; k < numberOfData; k++) {
            for (int j = k + 1; j < numberOfData; j++) {
                if (startDistances[k] == startDistances[j]) {
                    c++;
                    break;
                }
            }
        }
        if (c == 0) {
            printf("Found it pivot must be %d\n", i);
            break;
        }
        c = 0;
        if (i == numberOfData - 1) {
            printf("Dataset must be changed\n");
        }
    }*/

    MPI_Init(&argc, &argv);

    int world_size;
    int world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    //Splitting the datas to processes

    double dataPortion[numberOfData / world_size][numberOfCol];

    for (int i = 0; i < numberOfData / world_size; i++) {
        for (int j = 0; j < numberOfCol; j++) {
            dataPortion[i][j] = data[world_rank * (numberOfData / world_size) + i][j];
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 0; i < numberOfData / world_size; i++) {
        for (int j = 0; j < numberOfCol; j++) {
            dataPortion[i][j] = data[world_rank * (numberOfData / world_size) + i][j];
        }
    }

    //Call of the distributeByMedian function

    double* pivot;
    pivot = (double*)malloc(numberOfCol * sizeof(double));

    MPI_Barrier(MPI_COMM_WORLD);

    distributeByMedian(0, MPI_COMM_WORLD, numberOfData, numberOfCol, numberOfRows, dataPortion, pivot);

    MPI_Barrier(MPI_COMM_WORLD);

    //Autocheck

    double finalDistances = 0;
    double minDistance = 0;
    double maxDistance = 0;

    for (int i = 0; i < numberOfData / world_size; i++) {
        finalDistances = euclidianDistance(pivot, dataPortion[i], numberOfCol);
        if ((finalDistances < minDistance) || (i == 0)) {
            minDistance = finalDistances;
        }
        if (finalDistances > maxDistance || (i == 0)) {
            maxDistance = finalDistances;
        }
    }
    double* minDistancesArr = NULL;
    double* maxDistancesArr = NULL;

    if (world_rank == 0) {
        minDistancesArr = (double*)malloc(world_size * sizeof(double));
        maxDistancesArr = (double*)malloc(world_size * sizeof(double));
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(&minDistance, 1, MPI_DOUBLE, minDistancesArr, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&maxDistance, 1, MPI_DOUBLE, maxDistancesArr, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    if (world_rank == 0) {

        int flag;
        for (int i = 0; i < (world_size - 1); i++) {
            if (maxDistancesArr[i] > minDistancesArr[i + 1]) {
                printf("ERROR, rank %d\n max of %d is %.3f \n min of %d is %.3f\n", i, i, maxDistancesArr[i], i + 1, minDistancesArr[i + 1]);
                flag = 1;
            }
        }

        if (flag == 0) {
            printf("--------------------------\n\nCode Executed Successfully\n\n--------------------------\n");
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();

    // Calculate time it took
    if (world_rank == 0) {
        clock_gettime(CLOCK_REALTIME, &requestEnd);

        double accum = (requestEnd.tv_sec - requestStart.tv_sec)
            + (requestEnd.tv_nsec - requestStart.tv_nsec)
                / BILLION;
        printf("Time it took :%lf seconds\n", accum);
    }

    return 0;
}