#ifndef INSERT_NODE_THREADS
#define INSERT_NODE_THREADS

#include <stdio.h>
#include <stdlib.h>

#include "../../seq_functions/linked_list/insert_node.h"

void vp_tree_threads(struct node* root, struct node** nodes, int* node_counter, int d, int num_threads);

#endif
