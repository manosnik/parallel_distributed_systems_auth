#include "eucDist.h"
#include <math.h>
#include <stdio.h>

double eucDist(double* center, double* point, int dimensions)
{

    double distance = 0;
    for (int i = 0; i < dimensions; i++) {
        distance += pow(point[i] - center[i], 2);
    }
    return distance;
}
