#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "../../seq_functions/linked_list/insert_node.h"
#include "inner_outer_threads.h"

typedef struct Struct {
    int id;
    int num_threads;
    struct node* root;
    int dim;
    double* distances;
    int* data_parts_size_inner;
    int* data_parts_size_outer;
    pthread_mutex_t* mutex;
} makeStruct;

typedef struct Struct3 {
    int id;
    int num_threads;
    struct node* root;
    int dim;
    double* distances;
    int* data_parts_size_inner;
    int* data_parts_size_outer;
    double** data_inner;
    double** data_outer;
    int* inner;
    int* outer;
    bool* flag;
    pthread_mutex_t* mtxin;
    pthread_mutex_t* mtxout;
} makeStruct3;

void* inner_outer_threads(void* args)
{
    makeStruct* arg = (makeStruct*)args;
    int range = arg->root->data_size / arg->num_threads;
    int start = arg->id * range;
    int end = start + range;
    if (arg->id == arg->num_threads - 1) {
        end = arg->root->data_size;
    }
    for (int i = start; i < end; i++) {
        if (arg->distances[i] <= arg->root->median_distance) {
            pthread_mutex_lock(arg->mutex);
            (*arg->data_parts_size_inner)++;
            (*arg->data_parts_size_outer)--;
            pthread_mutex_unlock(arg->mutex);
        }
    }
    /*
    // prepare conditions to wait for other threads

    pthread_mutex_lock(arg->lock);
    (*arg->lock_counter)++;
    // printf("BFOR: id: %d____lock counter is: %d\n", arg->id, *arg->lock_counter);
    pthread_mutex_unlock(arg->lock);

    // wait for all threads to finish

    pthread_mutex_lock(arg->lock);
    while (*arg->lock_counter < arg->num_threads) {
        // printf("WAIT: id: %d____lock counter is: %d\n", arg->id, *arg->lock_counter);
        pthread_cond_wait(arg->cond1, arg->lock);
    }
    pthread_mutex_unlock(arg->lock);

    pthread_mutex_lock(arg->lock);
    if ((*arg->lock_counter) == arg->num_threads) {
        // printf("FREE: id: %d____lock counter is: %d\n", arg->id, *arg->lock_counter);
        *arg->lock_counter = *arg->lock_counter + 1;
        pthread_cond_broadcast(arg->cond1);
    }
    pthread_mutex_unlock(arg->lock);
    */
    return NULL;
}

void* inner_outer_split(void* args)
{
    makeStruct3* arg = (makeStruct3*)args;
    int range = arg->root->data_size / arg->num_threads;
    int start = arg->id * range;
    int end = start + range;
    if (arg->id == arg->num_threads - 1) {
        end = arg->root->data_size;
    }
    // assign each data point to inner or outer node
    for (int i = start; i < end; i++) {
        if (arg->distances[i] < arg->root->median_distance) {
            pthread_mutex_lock(arg->mtxin);
            for (int j = 0; j < arg->dim; j++) {
                arg->data_inner[*(arg->inner)][j] = arg->root->data[i][j];
            }
            (*arg->inner) = (*arg->inner) + 1;
            pthread_mutex_unlock(arg->mtxin);

        } else if (arg->distances[i] > arg->root->median_distance) {
            pthread_mutex_lock(arg->mtxout);
            for (int j = 0; j < arg->dim; j++) {
                arg->data_outer[*(arg->outer)][j] = arg->root->data[i][j];
            }
            (*arg->outer) = (*arg->outer) + 1;
            pthread_mutex_unlock(arg->mtxout);
        } else if ((arg->distances[i] > arg->root->median_distance) && (*(arg->flag) == 0)) {
            pthread_mutex_lock(arg->mtxin);
            for (int j = 0; j < arg->dim; j++) {
                arg->data_inner[*(arg->inner)][j] = arg->root->data[i][j];
            }
            (*arg->inner) = (*arg->inner) + 1;
            *arg->flag = 1;
            pthread_mutex_unlock(arg->mtxin);
        } else if ((arg->distances[i] > arg->root->median_distance) && (*(arg->flag) == 1)) {
            pthread_mutex_lock(arg->mtxout);
            for (int j = 0; j < arg->dim; j++) {
                arg->data_outer[*(arg->outer)][j] = arg->root->data[i][j];
            }
            (*arg->outer) = (*arg->outer) + 1;
            *arg->flag = 0;
            pthread_mutex_unlock(arg->mtxout);
        }
    }
    return NULL;
}