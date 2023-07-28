#include "insert_node.h"
#include "../../helpers/eucDist.h"
#include "../median.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

struct node* create_node(double** data, int id_vp, double* vp, int size, int d)
{
    struct node* new_node = (struct node*)malloc(sizeof(struct node));

    new_node->vp = (double*)malloc(d * sizeof(double));
    new_node->data = (double**)malloc(size * sizeof(double*));
    for (int i = 0; i < size; i++) {
        new_node->data[i] = (double*)malloc(d * sizeof(double));
    }
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < d; j++) {
            new_node->data[i][j] = data[i][j];
        }
    }
    // new_node->data = data;

    new_node->id = id_vp;

    new_node->vp = vp;

    new_node->data_size = size;

    new_node->inner = NULL;
    new_node->outer = NULL;

    new_node->prev = NULL;
    return new_node;
}

void vp_tree(struct node* root, struct node** nodes, int* node_counter, int d)
{
    if (root->data_size <= 1) {
        return;
    }

    // finding the median distance

    double* distances = malloc(root->data_size * sizeof(double));
    double* distances_copy = malloc(root->data_size * sizeof(double));

    for (int i = 0; i < root->data_size; i++) {
        distances[i] = eucDist(root->vp, root->data[i], d);
        distances_copy[i] = eucDist(root->vp, root->data[i], d);
    }

    root->median_distance = median(distances_copy, root->data_size, (root->data_size / 2) + 1);
    free(distances_copy);

    int data_parts_size_inner = 0;

    int data_parts_size_outer = 0;

    bool flag = 0;

    for (int i = 0; i < root->data_size; i++) {
        if (distances[i] < root->median_distance) {
            data_parts_size_inner++;
        } else if ((distances[i] == root->median_distance) && flag == 0) {
            data_parts_size_inner++;
            flag = 1;
        } else if ((distances[i] == root->median_distance) && flag == 1) {
            flag = 0;
        }
    }

    data_parts_size_outer = root->data_size - data_parts_size_inner;

    // assign each data point to inner or outer node

    double** data_inner = (double**)malloc(data_parts_size_inner * sizeof(double*));
    for (int i = 0; i < data_parts_size_inner; i++) {
        data_inner[i] = (double*)malloc(d * sizeof(double));
    }

    double** data_outer = (double**)malloc(data_parts_size_outer * sizeof(double*));
    for (int i = 0; i < data_parts_size_outer; i++) {
        data_outer[i] = (double*)malloc(d * sizeof(double));
    }

    int inner = 0, outer = 0;

    flag = 0;

    for (int i = 0; i < root->data_size; i++) {
        if (distances[i] < root->median_distance) {
            for (int j = 0; j < d; j++) {
                data_inner[inner][j] = root->data[i][j];
            }
            inner++;
        } else if ((distances[i] == root->median_distance) && flag == 0) {
            for (int j = 0; j < d; j++) {
                data_inner[inner][j] = root->data[i][j];
            }
            inner++;
            flag = 1;

        } else if ((distances[i] == root->median_distance) && flag == 1) {
            for (int j = 0; j < d; j++) {
                data_outer[outer][j] = root->data[i][j];
            }
            outer++;
            flag = 0;
        }

        else if (distances[i] > root->median_distance) {
            for (int j = 0; j < d; j++) {
                data_outer[outer][j] = root->data[i][j];
            }
            outer++;
        }
    }

    // create inner and outer nodes
    nodes[*node_counter] = create_node(data_inner, 0, data_inner[0], data_parts_size_inner, d);
    root->inner = nodes[*node_counter];
    nodes[*node_counter]->prev = root;
    (*node_counter)++;
    nodes[*node_counter] = create_node(data_outer, 0, data_outer[0], data_parts_size_outer, d);
    root->outer = nodes[*node_counter];
    nodes[*node_counter]->prev = root;
    (*node_counter)++;
}

void print_info(struct node* node, int d)
{
    printf("\n\n--------------------\n\n");
    printf("node id: %d\n", node->id);
    printf("node vp: ");
    for (int i = 0; i < d; i++) {
        printf("%f\t", node->vp[i]);
    }
    printf("\n");
    printf("node data size: %d\nmedian distance: %f\n\n-- -- -\tdata\t-- -- -\n",
        node->data_size,
        node->median_distance);
    for (int i = 0; i < node->data_size; i++) {
        for (int j = 0; j < d; j++) {
            printf("%f ", node->data[i][j]);
        }
        printf("\n");
    }
    if (node->inner != NULL) {
        printf("we got an inner child!\n");
    }
    if (node->outer != NULL) {
        printf("we got an outer child!\n");
    }
    printf("\n\n--------------------\n\n");
}