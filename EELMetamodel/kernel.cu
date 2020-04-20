/*
File: kernel.cu
Date: 2-24-2020
Author: Joe Haun
Purpose: This file is for the purposes of writing Metamodel test code. 
*/

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "EELSurrogate.cu"
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <complex>

void printMatrix(double inArray[], int dimension);

int main()
{
    const int arraySize = 4;
    double a[arraySize] = { 1, 2, 3, 4 };
    double b[arraySize] = { 10, 20, 30, 40};
    double c[arraySize] = { 5 };
    cudaError_t cudaStatus = cudaSuccess;
    
    // Add vectors in parallel.
    //cudaStatus = cudaSetup(c, a, b, sqrt(arraySize));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "Kernel failed!");
        return 1;
    }

    //Print first input matrix
    printMatrix(a, sqrt(arraySize));
    //Print second input matrix
    printMatrix(b, sqrt(arraySize));
    //Print output matrix
    printMatrix(c, sqrt(arraySize));

    // cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
    cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceReset failed!");
        return 1;
    }
   
    return 0;
}


void printMatrix(double inArray[], int dimension) {
    int index = 0;
    std::cout << "{";
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            index = i + (j * dimension);
            if (j + 1 < dimension) {
                std::cout << inArray[index] << ", ";
            }
            else {
                std::cout << inArray[index];
            }
        }
        if (i + 1 < dimension) {
            std::cout << "; ";
        }
        else {
            std::cout << "}" << std::endl;
        }
    }
    return;
}
