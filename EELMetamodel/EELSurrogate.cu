/*
File: EELSurrogate.cu
Date: 2-24-2020
Author: Joe Haun
Purpose: This file contains the EEL Surrogate class.  
*/

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <iostream>

//Function sets up and performs calculations on the GPU
cudaError_t cudaSetup(double* c, double* a, double* b, unsigned int size);

//Function calculates the covariance of the given matrices
void calculateCovariance(double outputMatrix[], double firstInput[], double secondInput[], int dimension);

//Function calculates the average of the given matrix
double calculateAverage(double inputMatrix[], int dimension);

//Function prints the given matrix on the GPU
__global__ void printWithGPU(double* matrix, int dimension);

//Function calculates the covariance of the given matrices and places them into a new output matrix
__global__ void covarianceHelper(double* outputMatrix, double* X, double meanX, double* Y, double meanY, int dimension);

//Function 
__global__ void nonDiagNormalize(double* firstInput, double* secondInput, int n, int i);
__global__ void diagNormalize(double* firstInput, double* secondInput, int n, int i);
__global__ void gaussJordan(double* firstInput, double* secondInput, int n, int i);
__global__ void setZero(double* input, int n, int i);

// Helper function for using CUDA to add vectors in parallel.
cudaError_t cudaSetup(double* c, double* a, double* b, unsigned int size)
{
    double* deviceFirstInput = 0;
    double* deviceSecondInput = 0;
    double* deviceOutput = 0;
    cudaError_t cudaStatus;
    int memoryAllocationSize = pow(size, 2);

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }

    // Allocate GPU buffers for three vectors (two input, one output)    .
    cudaStatus = cudaMalloc((void**)&deviceOutput, memoryAllocationSize * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&deviceFirstInput, memoryAllocationSize * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&deviceSecondInput, memoryAllocationSize * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    // Copy input vectors from host memory to GPU buffers.
    cudaStatus = cudaMemcpy(deviceFirstInput, a, memoryAllocationSize * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(deviceSecondInput, b, memoryAllocationSize * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    // Launch a kernel on the GPU with one thread for each element.
    //addKernel << <size, size >> > (deviceOutput, deviceFirstInput, deviceSecondInput);
    double firstInputMean = calculateAverage(a, size);
    double secondInputMean = calculateAverage(b, size);
    covarianceHelper <<<size, size>>> (deviceOutput, deviceFirstInput, firstInputMean, deviceSecondInput, secondInputMean, size);
    //printWithGPU << <size, size >> > (deviceFirstInput, size);
    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }

    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n\n", cudaStatus);
        fprintf(stderr, "MSG: %s\n\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }

    // Copy output vector from GPU buffer to host memory.
    cudaStatus = cudaMemcpy(c, deviceOutput, size * size * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

Error:
    cudaFree(deviceOutput);
    cudaFree(deviceFirstInput);
    cudaFree(deviceSecondInput);

    return cudaStatus;
}

__global__ void printWithGPU(double* matrix, int dimension) {

    printf("Hello from <%d, %d>: %f\n", blockIdx.x, threadIdx.x, matrix[blockIdx.x + (threadIdx.x * dimension)]);
    return;
}

__global__ void covarianceHelper(double* outputMatrix, double* X, double meanX, double* Y, double meanY, int dimension) {
    int i = threadIdx.x;
    int j = blockIdx.x;
    int index = i + (j * dimension);
    outputMatrix[index] = (X[index] - meanX) * (Y[index] - meanY) / dimension;
    return;
}

double calculateAverage(double inputMatrix[], int dimension) {
    //Calculates the average value of the given matrix. Likely a way to do this with CUDA.
    double meanVal;
    double rowSum = 0;
    double colSum = 0;
    for (int i = 0; i < dimension; i++) {
        colSum = 0;
        for (int j = 0; j < dimension; j++) {
            colSum += inputMatrix[i + (j * dimension)];
        }
        rowSum += colSum / dimension; //Divide column sum by dimension as there should be that many values present
    }

    meanVal = rowSum / dimension; //Divide row sum by dimension as there should be that many values present. 
    //This achieves a division by the total number of values AKA dimension*dimension.
    return meanVal;
}

__global__ void nonDiagNormalize(double* firstInput, double* secondInput, int n, int i) { 

    return; 
}