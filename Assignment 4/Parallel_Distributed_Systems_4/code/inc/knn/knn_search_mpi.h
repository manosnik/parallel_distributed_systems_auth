#ifndef KNN_SEARCH_MPI
#define KNN_SEARCH_MPI

#include "../seq_functions/linked_list/insert_node.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

void knn_search_mpi(int chief_rank, struct node** nodes, int N, int k, int d);

#endif