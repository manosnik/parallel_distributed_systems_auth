#include <math.h>
#include <stdio.h>

int tree_info_calc(int N)
{
    int height = 0;
    int nodes = 0;

    // height is upper log2(N) with math.h log2(N)
    height = (int)ceil(log2(N));

    //  node is sum of 2^height
    for (int i = 0; i <= height; i++) {
        nodes += pow(2, i);
    }

    return nodes;
}