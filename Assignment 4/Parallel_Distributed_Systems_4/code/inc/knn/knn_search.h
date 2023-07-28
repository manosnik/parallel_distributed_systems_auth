#ifndef KNN_SEARCH
#define KNN_SEARCH

#include "../seq_functions/linked_list/insert_node.h"

#include <stdio.h>
#include <stdlib.h>

void knn_search(struct node** nodes, int dim, double** knn, int k, int knn_pos);

#endif