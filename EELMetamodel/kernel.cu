/*
File: kernel.cu
Date: 2-24-2020
Modified: 8-2-2020
Author: Joe Haun
Purpose: This file is for the purposes of writing Metamodel test code. 
*/

#include "EELSurrogate.cu"

int main(int argc, char* argv)
{
    const int arraySize = 9;
    double dummyDesignSite[arraySize] = { 2, 1, 3, 5, 2, 3, 5, 2, 9};
    double dummyTestSite[arraySize] = { 10, 11, 0, 0};
    double dummyDesignSiteValues[arraySize] = { 1, 3, 4, 22, 9, 0, 8, 2, 11 };
    double result = -1;

    double theta = 1;
    double aVal = 1;
    double variance = 0;
    cudaError_t cudaStatus = cudaSuccess;
    
    printf("MAIN\n\nDesign Site\n");
    printMatrix(dummyDesignSite, sqrt(arraySize));
    printf("\nDesign Site Values\n");
    printMatrix(dummyDesignSiteValues, sqrt(arraySize));

    //Perform metamodel Kriging operations
    result = metamodelSetup(sqrt(arraySize), theta, variance, aVal, dummyDesignSite, dummyTestSite, dummyDesignSiteValues);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "Kernel failed!");
        return 1;
    }
    /*
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
    */
    // cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
    cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceReset failed!");
        return 1;
    }
   
    return 0;
}



