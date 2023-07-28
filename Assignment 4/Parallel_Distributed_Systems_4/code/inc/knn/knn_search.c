#include "knn_search.h"
#include "../helpers/eucDist.h"

#include "../seq_functions/median.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void knn_search(struct node** nodes, int dim, double** knn, int k, int knn_pos)
{
    struct node* current_node = nodes[0];

    double distance = 0;

    // find the leaf node

    while (current_node->inner != NULL) {
        if (k + 1 > current_node->data_size) {
            current_node = current_node->prev;
            break;
        }
        distance = eucDist(current_node->vp, nodes[0]->data[knn_pos], dim);
        if (distance < current_node->median_distance) {
            current_node = current_node->inner;
        } else {
            current_node = current_node->outer;
        }
    }

    // find the k nearest neighbors
    double temp_distance;
    double* temp_distances_arr = (double*)malloc(k * sizeof(double));

    for (int i = 0; i < current_node->data_size; i++) {
        if (i < k) {
            for (int j = 0; j < dim; j++) {
                knn[i][j] = current_node->data[i][j];
            }
            if (eucDist(current_node->data[i], nodes[0]->data[knn_pos], dim) == 0) {
                temp_distances_arr[i] = INFINITY;
            } else {
                temp_distances_arr[i] = eucDist(current_node->data[i], nodes[0]->data[knn_pos], dim);
            }

        } else {
            temp_distance = eucDist(current_node->data[i], nodes[0]->data[knn_pos], dim);
            for (int j = 0; j < k; j++) {

                if ((temp_distance < temp_distances_arr[j]) && (temp_distance != 0.0)) {
                    temp_distances_arr[j] = temp_distance;
                    for (int k = 0; k < dim; k++) {
                        knn[j][k] = current_node->data[i][k];
                    }
                    break;
                }
            }
        }
    }

    double max_distance = 0;
    double vp_distance = 0;
    double rad_distance = 0;

    while (current_node->data_size <= nodes[0]->data_size) {
        // find the max distance of the k nearest neighbors
        max_distance = temp_distances_arr[find_max(temp_distances_arr, k)];

        vp_distance = eucDist(current_node->vp, nodes[0]->data[knn_pos], dim);

        rad_distance = current_node->median_distance - vp_distance;

        // compare the max distance of the k nearest neighbors with the radius of the current node

        if (rad_distance < max_distance) {
            // if the max distance is smaller than the radius, then we have to search the other side of the tree
            // find the other side of the tree
            if (current_node->prev->inner == current_node) {
                current_node = current_node->prev->outer;
            } else {
                current_node = current_node->prev->inner;
            }
            // find the k nearest neighbors
            for (int i = 0; i < current_node->data_size; i++) {

                temp_distance = eucDist(current_node->data[i], nodes[0]->data[knn_pos], dim);
                for (int j = 0; j < k; j++) {
                    if (temp_distance < temp_distances_arr[j] && temp_distance != 0.0) {
                        temp_distances_arr[j] = temp_distance;
                        for (int p = 0; p < dim; p++) {
                            knn[j][p] = current_node->data[i][p];
                        }
                        break;
                    }
                }
            }
            current_node = current_node->prev;
        } else {
            break;
        }
        if (current_node->data_size == nodes[0]->data_size) {
            break;
        }
    }
    if (knn_pos == 7) {
        for (int i = 0; i < k; i++) {
            if (temp_distances_arr[i] == 0) {
            }
            printf("%f\t", temp_distances_arr[i]);
        }
    }

    return;
}