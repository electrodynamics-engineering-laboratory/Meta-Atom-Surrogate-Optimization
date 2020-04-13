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
cudaError_t cudaSetup(double* c, double* a, double* b, double* identityMatrix, unsigned int size);

//Function calculates the average of the given matrix
double calculateAverage(double inputMatrix[], int dimension);

//Function prints the given matrix on the GPU
__global__ void printWithGPU(double* matrix, int dimension);

//Function calculates the covariance of the given matrices and places them into a new output matrix
__global__ void calculateCovariance(double* outputMatrix, double* X, double meanX, double* Y, double meanY, int dimension);

//Functions are used to calculate the inverse of a given matrix using Gauss-Jordan elimination method
__global__ void nonDiagNormalize(double* inputMatrix, double* identityMatrix, int n, int i);
__global__ void diagNormalize(double* inputMatrix, double* identityMatrix, int n, int i);
__global__ void gaussJordan(double* inputMatrix, double* identityMatrix, int n, int i);
__global__ void setZero(double* inputMatrix, double* identityMatrix, int n, int i);

//Function calculates the Kriging Estimator for the Global Mean
__global__ void calculateKrigingEstimator(double* result, double* covarianceMatrix, double* sampleMatrix);

//Functions calculate the dot product of two vectors, overloaded in cases where only a scalar is returned
__global__ void calculateDotProduct(double* result, double* inputOne, int xOne, int yOne, double* inputTwo, int xTwo, int yTwo);
__global__ void calculateDotProduct(double result, double* inputOne, int xOne, int yOne, double* inputTwo, int xTwo, int yTwo);

// Helper function for using CUDA to add vectors in parallel.
cudaError_t cudaSetup(double* c, double* a, double* b, double* identityMatrix, unsigned int size)
{
    double* deviceFirstInput = 0;
    double* deviceSecondInput = 0;
    double* deviceOutput = 0;
    double* deviceIdentityMatrix = 0;
    cudaError_t cudaStatus;
    int memoryAllocationSize = pow(size, 2);

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }

    // Allocate GPU buffers for three vectors (two input, one output)    .
    cudaStatus = cudaMalloc((void**)&deviceIdentityMatrix, memoryAllocationSize * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

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

    cudaStatus = cudaMemcpy(deviceIdentityMatrix, identityMatrix, memoryAllocationSize * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    // Launch a kernel on the GPU with one thread for each element.
    double firstInputMean = calculateAverage(a, size);
    double secondInputMean = calculateAverage(b, size);
    calculateCovariance <<<size, size>>> (deviceOutput, deviceFirstInput, firstInputMean, deviceSecondInput, secondInputMean, size);
    
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
    cudaStatus = cudaMemcpy(c, deviceOutput, memoryAllocationSize * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

Error:
    cudaFree(deviceOutput);
    cudaFree(deviceFirstInput);
    cudaFree(deviceSecondInput);
    cudaFree(deviceIdentityMatrix);

    return cudaStatus;
}

__global__ void printWithGPU(double* matrix, int dimension) {

    printf("Hello from <%d, %d>: %f\n", blockIdx.x, threadIdx.x, matrix[blockIdx.x + (threadIdx.x * dimension)]);
    return;
}

__global__ void calculateCovariance(double* outputMatrix, double* X, double meanX, double* Y, double meanY, int dimension) {
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

__global__ void nonDiagNormalize(double* inputMatrix, double* identityMatrix, int n, int i) { 
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if ((x < n) && (y < n)) {
        if ((x == i) && (x != y)) {
            identityMatrix[x * n + y] /= inputMatrix[i * n + i];
            inputMatrix[x * n + y] /= inputMatrix[i * n + i];
        }
    }
    return; 
}

__global__ void diagNormalize(double* inputMatrix, double* identityMatrix, int n, int i) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    
    if ((x < n) && (y < n)) {
        if ((x == i) && (x != y)) {
            identityMatrix[x * n + y] /= inputMatrix[i * n + i];
            inputMatrix[x * n + y] /= A[i * n + i];
        }
    }

    return;
}

__global__ void gaussJordan(double* inputMatrix, double* identityMatrix, int n, int i) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if ((x < n) && (y < n)) {
        if (x != i) {
            identityMatrix[x * n + y] -= identityMatrix[i * n + y] * inputMatrix[x * n + i];
            if (y != i) {
                inputMatrix[x * n + y] -= A[i * n + y] * inputMatrix[x * n + i];
            }
        }
    }

    return;
}

__global__ void setZero(double* inputMatrix, double* identityMatrix, int n, int i) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if ((x < n) && (y < n)) {
        if (x != i) {
            if (y != i) {
                inputMatrix[x * n + y] = 0;
            }
        }
    }

    return;
}

__global__ void calculateKrigingEstimator(double* result, double* covarianceMatrix, double* sampleMatrix) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;



    return;
}

