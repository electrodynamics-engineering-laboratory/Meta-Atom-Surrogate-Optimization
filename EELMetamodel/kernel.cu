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


int main()
{
    const int arraySize = 4;
    double dummyDesignSite[arraySize] = { 1, 2, 3, 4 };
    double dummyTestSite[arraySize] = { 10, 11, 0, 0};
    double dummyDesignSiteValues[arraySize] = { 6, 7, 8, 9 };
    double result = -1;

    double theta = 1;
    double aVal = 1;
    double variance = 0;
    cudaError_t cudaStatus = cudaSuccess;

    
    //Perform metamodel Kriging operations
    result = metamodelSetup(sqrt(arraySize), theta, variance, aVal, dummyDesignSite, dummyTestSite, dummyDesignSiteValues);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "Kernel failed!");
        return 1;
    }

    std::cout << "MAIN: Matrices after calculations." << std::endl;
    //Print first input matrix
    printf("MAIN: designSite\n");
    printMatrix(dummyDesignSite, sqrt(arraySize));
    //Print second input matrix
    printf("MAIN: testSite\n");
    printMatrix(dummyTestSite, sqrt(arraySize));
    //Print output matrix
    printf("MAIN: designSiteValues\n");
    printMatrix(dummyDesignSiteValues, sqrt(arraySize));

    // cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
    cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceReset failed!");
        return 1;
    }
   
    return 0;
}



