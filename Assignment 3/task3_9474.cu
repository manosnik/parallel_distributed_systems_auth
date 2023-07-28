/*-----------------------------------------------------------
The program has been run at Google Colab with the following orders
-------------------------------------------------------------
!apt-get --purge remove cuda nvidia* libnvidia-*
!dpkg -l | grep cuda- | awk '{print $2}' | xargs -n1 dpkg --purge
!apt-get remove cuda-*
!apt autoremove
!apt-get update
--------------------------------------------------------------
!wget https://developer.nvidia.com/compute/cuda/9.2/Prod/local_installers/cuda-repo-ubuntu1604-9-2-local_9.2.88-1_amd64 -O cuda-repo-ubuntu1604-9-2-local_9.2.88-1_amd64.deb
!dpkg -i cuda-repo-ubuntu1604-9-2-local_9.2.88-1_amd64.deb
!apt-key add /var/cuda-repo-9-2-local/7fa2af80.pub
!apt-get update
!apt-get install cuda-9.2
--------------------------------------------------------------
!nvcc --version
--------------------------------------------------------------
%load_ext nvcc_plugin
--------------------------------------------------------------
and at the AUTH's cluster with the known procedure*/

//%%cu necessary in Google Colab

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define BILLION 1E9

#define N 2000 //N=M
#define M 2000
#define K 31     //must be K%2=1
#define b 20     //max b around 10 because of the size of shared memory
#define BV1 4000 //BV1 should not be more than 1000x smaller than N*M
#define BV2 2000
#define BV3 2000 //1.The number of threads of every block in version 3 multiplied with b should not exceed the size of shared memory and \
                   //2.The number of columns(=threadsPerBlock*b) that everuy block computes should follow .. M%Block Col==0

//function of serial code

int sign(int a, int t, int c, int d, int e)
{

    int temp = 0;
    temp = a + t + c + d + e;
    if (temp > 0)
    {
        return 1;
    }
    else if (temp < 0)
    {
        return -1;
    }
    else if (temp == 0)
    {
        printf("Error in sign function\n");
        return 0;
    }
    return 7;
}

//function of v1

__global__ void isingModel_v1(int *initialArray, int *secondaryArray)
{
    int index = threadIdx.x + blockIdx.x * blockDim.x;

    int i = index / M;
    int j = index % N;

    secondaryArray[i * N + j] = initialArray[((i - 1 + N) % N) * N + (j % M)] + initialArray[(i % N) * N + ((j - 1 + M) % M)] + initialArray[(i % N) * N + (j % M)] + initialArray[((i + 1) % N) * N + (j % M)] + initialArray[(i % N) * N + ((j + 1) % M)];

    if (secondaryArray[i * N + j] > 0)
    {
        secondaryArray[i * N + j] = 1;
    }
    else if (secondaryArray[i * N + j] < 0)
    {
        secondaryArray[i * N + j] = -1;
    }
    else if (secondaryArray[i * N + j] == 0)
    {
        secondaryArray[i * N + j] = 3;
    }

    __syncthreads();
}

//function v2

__global__ void isingModel_v2(int *initialArray, int *secondaryArray)
{
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    int rowStart = (int)(index / (int)(M / b)) * b;
    int colStart = (int)(index % (int)(M / b)) * b;

    for (int i = rowStart; i < rowStart + b & i < N; i++)
    {
        for (int j = colStart; j < colStart + b & j < M; j++)
        {
            secondaryArray[i * N + j] = initialArray[((i - 1 + N) % N) * N + (j % M)] + initialArray[(i % N) * N + ((j - 1 + M) % M)] + initialArray[(i % N) * N + (j % M)] + initialArray[((i + 1) % N) * N + (j % M)] + initialArray[(i % N) * N + ((j + 1) % M)];

            if (secondaryArray[i * N + j] > 0)
            {
                secondaryArray[i * N + j] = 1;
            }
            else if (secondaryArray[i * N + j] < 0)
            {
                secondaryArray[i * N + j] = -1;
            }
            else if (secondaryArray[i * N + j] == 0)
            {
                secondaryArray[i * N + j] = 3;
            }
        }
    }
    __syncthreads();
}

//function of v3

__global__ void isingModel_v3(int *initialArray, int *secondaryArray)
{
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    int rowStart = (int)(index / (int)(N / b)) * b;
    int colStart = (int)(index % (int)(M / b)) * b;

    int threadsPerBlock = ((N * M) / (b * b)) / BV3;
    int cudaArrJ = threadsPerBlock * b;
    int cudaArrI = M / (threadsPerBlock * b);

    __shared__ int cudaSharedArr[M * N / BV3];

    for (int i = rowStart; i < rowStart + b & i < N; i++)
    {
        for (int j = colStart; j < colStart + b & j < M; j++)
        {

            cudaSharedArr[(i % b) * cudaArrJ + j % cudaArrJ] = initialArray[i * M + j];
        }
    }

    __syncthreads();

    int a, c, d, e;
    for (int i = rowStart; i < rowStart + b & i < N; i++)
    {
        for (int j = colStart; j < colStart + b & i < M; j++)
        {

            if (i - 1 == rowStart - 1)
            {
                a = initialArray[((i - 1 + N) % N) * M + (j % M)];
            }
            else
            {
                a = cudaSharedArr[(((i - 1 + N) % N) % b) * cudaArrJ + j % cudaArrJ];
            }

            if (i + 1 == rowStart + b)
            {
                e = initialArray[((i + 1) % N) * M + (j % M)];
            }
            else
            {
                e = cudaSharedArr[(((i + 1) % N) % b) * cudaArrJ + j % cudaArrJ];
            }

            if (threadIdx.x == 0 && j == colStart)
            {
                c = initialArray[(i % N) * M + ((j - 1 + M) % M)];
            }
            else
            {
                c = cudaSharedArr[(i % b) * cudaArrJ + ((j - 1 + M) % M) % cudaArrJ];
            }

            if (threadIdx.x == blockDim.x - 1 && j == colStart + b - 1)
            {

                d = initialArray[(i % N) * M + ((j + 1) % M)];
            }
            else
            {
                d = cudaSharedArr[(i % b) * cudaArrJ + ((j + 1) % M) % cudaArrJ];
            }

            secondaryArray[i * M + j] = a + c + cudaSharedArr[(i % b) * cudaArrJ + j % cudaArrJ] + e + d;

            if (secondaryArray[i * M + j] > 0)
            {
                secondaryArray[i * M + j] = 1;
            }
            else if (secondaryArray[i * M + j] < 0)
            {
                secondaryArray[i * M + j] = -1;
            }
            else if (secondaryArray[i * M + j] == 0)
            {
                secondaryArray[i * M + j] = 3;
            }
        }
    }
}

int main(int argc, char *argv[])
{

    const int size = N * M * sizeof(int);

    // declare the initial and the secondry 2D array
    int **initialArray_v0 = (int **)malloc(N * sizeof(int *));
    int **secondaryArray_v0 = (int **)malloc(N * sizeof(int *));

    for (int i = 0; i < N; i++)
    {
        initialArray_v0[i] = (int *)malloc(M * sizeof(int));
        secondaryArray_v0[i] = (int *)malloc(M * sizeof(int));
    }

    // initialize the fist array for all versions
    srand(time(0));
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            initialArray_v0[i][j] = rand() % 2;
            if (initialArray_v0[i][j] == 0)
            {
                initialArray_v0[i][j] = -1;
            }
        }
    }

    int *initialArray_v1 = (int *)malloc(N * M * sizeof(int *));
    int *initialArray_v2 = (int *)malloc(N * M * sizeof(int *));
    int *initialArray_v3 = (int *)malloc(N * M * sizeof(int *));

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            initialArray_v1[i * M + j] = initialArray_v0[i][j];
            initialArray_v2[i * M + j] = initialArray_v0[i][j];
            initialArray_v3[i * M + j] = initialArray_v0[i][j];
        }
    }
    //Checking that the 4 versions begin from the same starting point
    int fd = 0;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            if (initialArray_v0[i][j] != initialArray_v1[i * N + j] || initialArray_v0[i][j] != initialArray_v2[i * N + j] || initialArray_v0[i][j] != initialArray_v3[i * N + j])
            {
                printf("Erorr\n");
                fd = 1;
                break;
            }
        }
    }
    if (fd == 0)
    {
        printf("4 identical matrices\n");
    }

    //Start the clock

    struct timespec requestStart, requestEnd;
    clock_gettime(CLOCK_REALTIME, &requestStart);

    // Ising model of serial
    for (int k = 0; k < K; k++)
    {

        for (int i = 0; i < N; i++)
        {

            for (int j = 0; j < M; j++)
            {
                secondaryArray_v0[i][j] = sign(initialArray_v0[(i - 1 + N) % N][j], initialArray_v0[i][(j - 1 + M) % M], initialArray_v0[i][j], initialArray_v0[(i + 1) % N][j], initialArray_v0[i][(j + 1) % M]);
            }
        }

        // Copying the pointers
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {
                initialArray_v0[i][j] = secondaryArray_v0[i][j];
            }
        }
    }

    // Calculate time it took

    clock_gettime(CLOCK_REALTIME, &requestEnd);

    double accum = (requestEnd.tv_sec - requestStart.tv_sec) + (requestEnd.tv_nsec - requestStart.tv_nsec) / BILLION;
    printf("Time it took v0 :%lf seconds\n", accum);

    //version 1

    int threads_v1 = N * M / BV1;

    int *cudaInitialArray_v1, *cudaSecondaryArray_v1;

    cudaMalloc((void **)&cudaInitialArray_v1, size);
    cudaMalloc((void **)&cudaSecondaryArray_v1, size);

    clock_gettime(CLOCK_REALTIME, &requestStart);
    cudaMemcpy(cudaInitialArray_v1, initialArray_v1, size, cudaMemcpyHostToDevice);
    for (int k = 0; k < K; k++)
    {
        if (k % 2 == 0)
        {
            isingModel_v1<<<BV1, threads_v1>>>(cudaInitialArray_v1, cudaSecondaryArray_v1);
        }
        else
        {
            isingModel_v1<<<BV1, threads_v1>>>(cudaSecondaryArray_v1, cudaInitialArray_v1);
        }
    }
    cudaMemcpy(initialArray_v1, cudaSecondaryArray_v1, size, cudaMemcpyDeviceToHost);
    // Calculate time it took

    clock_gettime(CLOCK_REALTIME, &requestEnd);

    accum = (requestEnd.tv_sec - requestStart.tv_sec) + (requestEnd.tv_nsec - requestStart.tv_nsec) / BILLION;
    printf("Time it took v1 :%lf seconds\n", accum);

    cudaFree(cudaInitialArray_v1);
    cudaFree(cudaSecondaryArray_v1);

    //version 2

    int threads_v2 = ((N * M) / (b * b)) / BV2;

    int *cudaInitialArray_v2, *cudaSecondaryArray_v2;

    cudaMalloc((void **)&cudaInitialArray_v2, size);
    cudaMalloc((void **)&cudaSecondaryArray_v2, size);

    clock_gettime(CLOCK_REALTIME, &requestStart);

    cudaMemcpy(cudaInitialArray_v2, initialArray_v2, size, cudaMemcpyHostToDevice);

    for (int k = 0; k < K; k++)
    {
        if (k % 2 == 0)
        {
            isingModel_v2<<<BV2, threads_v2>>>(cudaInitialArray_v2, cudaSecondaryArray_v2);
        }
        else
        {
            isingModel_v2<<<BV2, threads_v2>>>(cudaSecondaryArray_v2, cudaInitialArray_v2);
        }
    }

    cudaMemcpy(initialArray_v2, cudaSecondaryArray_v2, size, cudaMemcpyDeviceToHost);

    // Calculate time it took

    clock_gettime(CLOCK_REALTIME, &requestEnd);
    accum = (requestEnd.tv_sec - requestStart.tv_sec) + (requestEnd.tv_nsec - requestStart.tv_nsec) / BILLION;
    printf("Time it took v2 :%lf seconds\n", accum);

    cudaFree(cudaInitialArray_v2);
    cudaFree(cudaSecondaryArray_v2);

    //version 3

    int threads_v3 = ((N * M) / (b * b)) / BV3;

    int *cudaInitialArray_v3, *cudaSecondaryArray_v3;

    cudaMalloc((void **)&cudaInitialArray_v3, size);
    cudaMalloc((void **)&cudaSecondaryArray_v3, size);

    clock_gettime(CLOCK_REALTIME, &requestStart);
    cudaMemcpy(cudaInitialArray_v3, initialArray_v3, size, cudaMemcpyHostToDevice);

    for (int k = 0; k < K; k++)
    {
        if (k % 2 == 0)
        {

            isingModel_v3<<<BV3, threads_v3>>>(cudaInitialArray_v3, cudaSecondaryArray_v3);
        }
        else
        {
            isingModel_v3<<<BV3, threads_v3>>>(cudaSecondaryArray_v3, cudaInitialArray_v3);
        }
    }
    cudaMemcpy(initialArray_v3, cudaSecondaryArray_v3, size, cudaMemcpyDeviceToHost);
    // Calculate time it took

    clock_gettime(CLOCK_REALTIME, &requestEnd);
    accum = (requestEnd.tv_sec - requestStart.tv_sec) + (requestEnd.tv_nsec - requestStart.tv_nsec) / BILLION;
    printf("Time it took v3 :%lf seconds\n", accum);

    cudaFree(cudaInitialArray_v3);
    cudaFree(cudaSecondaryArray_v3);

    int fb = 0, errorsC = 0;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            if (initialArray_v0[i][j] != initialArray_v1[i * M + j] || initialArray_v0[i][j] != initialArray_v2[i * M + j] || initialArray_v0[i][j] != initialArray_v3[i * M + j])
            {
                if (errorsC < 10)
                {
                    fb = 1;
                    if (initialArray_v0[i][j] != initialArray_v1[i * M + j])
                    {
                        printf("Erorr in v1\n");
                    }
                    if (initialArray_v0[i][j] != initialArray_v2[i * M + j])
                    {
                        printf("Erorr in v2\n");
                    }
                    if (initialArray_v0[i][j] != initialArray_v3[i * M + j])
                    {
                        printf("Erorr in v3\n");
                    }
                    errorsC++;
                }
                break;
            }
        }
    }
    if (fb == 0)
    {
        printf("4 identical matrices at end\n");
    }

    return 0;
}