#ifndef MEDIAN
#define MEDIAN

#include <stdio.h>

void swap(double* a, double* b);
int partition(double* A, int left, int right);
double quickselect(double* A, int left, int right, int k);

double median(double* A, int n, int k);

int find_max(double* array, int n);

#endif