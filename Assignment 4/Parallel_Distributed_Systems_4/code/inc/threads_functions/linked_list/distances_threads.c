#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../helpers/eucDist.h"
#include "../../seq_functions/linked_list/insert_node.h"
#include "distances_threads.h"

typedef struct Struct {
    int id;
    int num_threads;
    struct node* root;
    int dim;
    double* distances;
    double* distances_copy;
} makeStruct;

void* distances_threads(void* args)
{
    makeStruct* arg = (makeStruct*)args;
    int range = arg->root->data_size / arg->num_threads;
    int start = arg->id * range;
    int end = start + range;
    if (arg->id == arg->num_threads - 1) {
        end = arg->root->data_size;
    }
    for (int i = start; i < end; i++) {
        arg->distances[i] = eucDist(arg->root->vp, arg->root->data[i], arg->dim);
        arg->distances_copy[i] = eucDist(arg->root->vp, arg->root->data[i], arg->dim);
    }
    return NULL;
}