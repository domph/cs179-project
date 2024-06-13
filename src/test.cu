#include "test.cuh"
#include <stdio.h>

__global__ void cuda_hello() {
    printf("Hello World from GPU!\n");
}

void call_cuda_hello() {
    cuda_hello<<<1, 1>>>();
}