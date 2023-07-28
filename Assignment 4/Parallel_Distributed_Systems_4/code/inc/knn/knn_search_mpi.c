#include "knn_search_mpi.h"
#include "../seq_functions/linked_list/insert_node.h"
#include "knn_search.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

void knn_search_mpi(int chief_rank, struct node** nodes, int N, int k, int d)
{
    int my_id = 0;
    int num_procs = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    int proc_data_length;

    if (my_id != num_procs - 1) {
        proc_data_length = (int)(N / num_procs);
    } else {
        proc_data_length = (int)(N / num_procs) + (N % num_procs);
    }

    double*** knn_proc = (double***)malloc(proc_data_length * sizeof(double**));
    for (int i = 0; i < proc_data_length; i++) {
        knn_proc[i] = (double**)malloc(k * sizeof(double*));
        for (int j = 0; j < k; j++) {
            knn_proc[i][j] = (double*)malloc(d * sizeof(double));
        }
    }
    int start_index = my_id * (int)(N / num_procs);

    for (int i = 0; i < proc_data_length; i++) {
        knn_search(nodes, d, knn_proc[i], k, i + start_index);
    }

    double* knn_proc_1D = (double*)malloc(proc_data_length * k * d * sizeof(double));
    // make 1D array
    for (int i = 0; i < proc_data_length; i++) {
        for (int j = 0; j < k; j++) {
            for (int l = 0; l < d; l++) {
                knn_proc_1D[i * k * d + j * d + l] = knn_proc[i][j][l];
            }
        }
    }

    double*** knn;
    double* knn_1D = NULL;
    int* receive_counts;
    int* displacements;

    if (my_id == chief_rank) {
        knn = (double***)malloc(N * sizeof(double**));
        for (int i = 0; i < N; i++) {
            knn[i] = (double**)malloc(k * sizeof(double*));
            for (int j = 0; j < k; j++) {
                knn[i][j] = (double*)malloc(d * sizeof(double));
            }
        }
        knn_1D = (double*)malloc(N * k * d * sizeof(double));

        receive_counts = (int*)malloc(num_procs * sizeof(int));
        displacements = (int*)malloc(num_procs * sizeof(int));

        for (int i = 0; i < num_procs; i++) {
            if (i != num_procs - 1) {
                receive_counts[i] = (int)(N / num_procs) * k * d;
            } else {
                receive_counts[i] = ((int)(N / num_procs) + (N % num_procs)) * k * d;
            }
            displacements[i] = i * (int)(N / num_procs) * k * d;
        }
    }

    MPI_Gatherv(knn_proc_1D, proc_data_length * k * d, MPI_DOUBLE, knn_1D, receive_counts, displacements, MPI_DOUBLE, chief_rank, MPI_COMM_WORLD);

    if (my_id == chief_rank) {
        // make 3D array
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < k; j++) {
                for (int l = 0; l < d; l++) {
                    knn[i][j][l] = knn_1D[i * k * d + j * d + l];
                }
            }
        }

        printf("\nknn data search for: ");
        for (int i = 0; i < d; i++) {
            printf("%lf ", nodes[0]->data[300][i]);
        }
        printf("\n");
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < d; j++) {
                printf("%lf ", knn[300][i][j]);
            }
            printf("\n");
        }
    }
    return;
}