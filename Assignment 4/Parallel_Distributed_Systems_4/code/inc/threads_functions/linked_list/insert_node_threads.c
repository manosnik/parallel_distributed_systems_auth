#include "insert_node_threads.h"
#include "../../helpers/eucDist.h"
#include "../../seq_functions/median.h"

#include "../../seq_functions/linked_list/insert_node.h"

#include "distances_threads.h"

#include "inner_outer_threads.h"

#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct Struct {
    int id;
    int num_threads;
    struct node* root;
    int dim;
    double* distances;
    double* distances_copy;
} makeStruct;

typedef struct Struct2 {
    int id;
    int num_threads;
    struct node* root;
    int dim;
    double* distances;
    int* data_parts_size_inner;
    int* data_parts_size_outer;
    pthread_mutex_t* mutex;
} makeStruct2;

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

void vp_tree_threads(struct node* root, struct node** nodes, int* node_counter, int d, int num_threads)
{
    // terminate condition

    if (root->data_size <= 1) {
        return;
    }

    makeStruct* arguments;
    pthread_t* threads;
    arguments = (makeStruct*)malloc(sizeof(makeStruct) * num_threads);

    threads = (pthread_t*)malloc(sizeof(pthread_t) * num_threads);

    // finding the median distance

    double* distances = malloc(root->data_size * sizeof(double));
    double* distances_copy = malloc(root->data_size * sizeof(double));

    for (int i = 0; i < num_threads; i++) {
        arguments[i].id = i;
        arguments[i].num_threads = num_threads;
        arguments[i].root = root;
        arguments[i].dim = d;
        arguments[i].distances = distances;
        arguments[i].distances_copy = distances_copy;
        pthread_create(&threads[i], NULL, distances_threads, (void*)&arguments[i]);
    }

    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    root->median_distance = median(distances_copy, root->data_size, (root->data_size / 2) + 1);
    free(distances_copy);
    int* data_parts_size_inner = malloc(sizeof(int));
    int* data_parts_size_outer = malloc(sizeof(int));

    *data_parts_size_inner = 0;
    *data_parts_size_outer = 0;

    bool* flag = malloc(sizeof(bool));
    *flag = 0;

    for (int i = 0; i < root->data_size; i++) {
        if (distances[i] < root->median_distance) {
            *data_parts_size_inner = *data_parts_size_inner + 1;
        } else if ((distances[i] == root->median_distance) && *flag == 0) {
            *data_parts_size_inner = *data_parts_size_inner + 1;
            *flag = 1;
        } else if ((distances[i] == root->median_distance) && *flag == 1) {
            *flag = 0;
        }
    }

    *data_parts_size_outer = root->data_size - *data_parts_size_inner;

    double** data_inner = (double**)malloc(sizeof(double*) * (*data_parts_size_inner));
    for (int i = 0; i < (*data_parts_size_inner); i++) {
        data_inner[i] = (double*)malloc(sizeof(double) * d);
    }
    double** data_outer = (double**)malloc(sizeof(double*) * (*data_parts_size_outer));
    for (int i = 0; i < (*data_parts_size_outer); i++) {
        data_outer[i] = (double*)malloc(sizeof(double) * d);
    }

    int* inner = malloc(sizeof(int));
    int* outer = malloc(sizeof(int));

    *inner = 0;
    *outer = 0;

    *flag = 0;

    makeStruct3* arguments3;
    pthread_t* threads3;

    arguments3 = (makeStruct3*)malloc(sizeof(makeStruct3) * num_threads);
    threads3 = (pthread_t*)malloc(sizeof(pthread_t) * num_threads);

    pthread_mutex_t* mtxin;
    mtxin = malloc(sizeof(pthread_mutex_t));
    pthread_mutex_init(mtxin, NULL);

    pthread_mutex_t* mtxout;
    mtxout = malloc(sizeof(pthread_mutex_t));
    pthread_mutex_init(mtxout, NULL);

    for (int i = 0; i < num_threads; i++) {
        arguments3[i].id = i;
        arguments3[i].num_threads = num_threads;
        arguments3[i].root = root;
        arguments3[i].dim = d;
        arguments3[i].distances = distances;
        arguments3[i].data_parts_size_inner = data_parts_size_inner;
        arguments3[i].data_parts_size_outer = data_parts_size_outer;
        arguments3[i].data_inner = data_inner;
        arguments3[i].data_outer = data_outer;
        arguments3[i].inner = inner;
        arguments3[i].outer = outer;
        arguments3[i].flag = flag;
        arguments3[i].mtxin = mtxin;
        arguments3[i].mtxout = mtxout;
        pthread_create(&threads3[i], NULL, inner_outer_split, (void*)&arguments3[i]);
    }
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads3[i], NULL);
    }
    pthread_mutex_destroy(mtxin);
    pthread_mutex_destroy(mtxout);

    // create inner and outer nodes
    nodes[*node_counter] = create_node(data_inner, 0, data_inner[0], *data_parts_size_inner, d);
    root->inner = nodes[*node_counter];
    nodes[*node_counter]->prev = root;
    (*node_counter)++;
    nodes[*node_counter] = create_node(data_outer, 0, data_outer[0], *data_parts_size_outer, d);
    root->outer = nodes[*node_counter];
    nodes[*node_counter]->prev = root;
    (*node_counter)++;
}
