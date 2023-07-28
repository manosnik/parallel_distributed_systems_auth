#ifndef INSERT_NODE
#define INSERT_NODE

#include <stdio.h>
#include <stdlib.h>

struct node {
    double* vp;
    int id;
    double median_distance;
    double** data;
    int data_size;
    struct node* inner;
    struct node* outer;
    struct node* prev;
};
struct node* create_node(double** data, int id_vp, double* vp, int size, int d);

void vp_tree(struct node* root, struct node** nodes, int* node_counter, int d);

void print_info(struct node* node, int d);

#endif
