/*
File: kernel.cu
Date: 2-24-2020
Modified: 1-31-2021
Author: Joe Haun
Purpose: This file is for the purposes of writing Metamodel test code. 
*/

#include <iostream>
#include "EELSurrogate.cu"

int main(int argc, char** argv)
{
    const int arraySize = 9;
    double dummyDesignSite[arraySize] = { 2, 1, 3, 5, 2, 3, 5, 2, 9};
    double dummyTestSite[arraySize] = { 10, 11, 0, 0};
    double dummyDesignSiteValues[arraySize] = { 1, 3, 4, 22, 9, 0, 8, 2, 11 };
    double result = -1;

    double theta = 0;
    double aVal = 0;
    double variance = 0;
    int dimension = 0;
    int headerLines = 0;
    cudaError_t cudaStatus = cudaSuccess;
    std::string filePath;

    printf("Enter a file with design site data: ");
    //std::cin >> filePath;
    filePath = "E:\\Documents\\PROJECTS\\SurrogateModeling\\Dogbone.csv";
    printf("Entered: %s\n", filePath.c_str());
    printf("Enter the dimensionality for the computation: ");
    std::cin >> dimension;
    if (dimension == 0) {
        printf("Please enter a valid integer dimension.\n");
        return 1;
    }
    printf("Entered: %d\n", dimension);
    printf("Enter the number of header lines in the file: ");
    std::cin >> headerLines;
    printf("Enter a value for Theta: ");
    std::cin >> theta;
    printf("Entered: %f\n", theta);
    printf("Enter a nugget effect value: ");
    std::cin >> aVal;
    printf("Entered: %f\n", aVal);
    printf("Enter a variance value: ");
    std::cin >> variance;
    printf("Entered: %f\n", variance);

    double* inputArray = readInputFile(filePath, dimension, headerLines);
    
    printf("Read Matrix\n");
    printMatrix(inputArray, 103, 10);
    //printMatrix(inputArray, 10);

    //Perform metamodel Kriging operations
    //result = metamodelSetup(sqrt(arraySize), theta, variance, aVal, dummyDesignSite, dummyTestSite, dummyDesignSiteValues);
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



