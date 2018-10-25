#include <iostream>
#include <fstream>
#include <omp.h>
#include <algorithm>
#include <stdio.h>

#define OUTPUT_FILE "out.txt"
#define MAX_ITER 20000

void calculate_parallel(int n, int num_t);