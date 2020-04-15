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
#include <complex>

//Begin CPU Function Definitions
/*
Name: metamodelSetup
Description: Function sets up and performs calculations on the GPU that relate to a metamodel
Input:  outputValue - (CURRENTLY UNKNOWN AND WILL CHANGE) a pointer to a value that will be used for output 
        dimension - an integer value representing the number of rows for all matrices and vectors and the number columns for all matrices
        theta - (THIS MIGHT BE AN INCORRECT DESCRIPTION) a value in radians that represents the angle between the two matrices, this value should only ever be real
        variance - the value of the variance between 
        designSite - a dimension-by-dimension matrix that represents the datapoints previously collected to be used to estimate a new value based on test site
        testSite - a one-by-dimension vector that represents the points at which a new value is to be calculated
        designSiteValues - a dimension-by-dimension matrix that represent the values at particular design sites. Each matrix has a 1-to-1 correlation of site to value (i.e. designSite[0] is paired with designSiteValues[0])
Output: a value that is any error that the GPU experienced upon attempting calculations
*/
//Function sets up and performs calculations on the GPU
cudaError_t metamodelSetup(std::complex<double>* outputValue, int dimension, std::complex<double> theta, std::complex<double> variance, std::complex<double> a, std::complex<double>* designSite, std::complex<double>* testSite, std::complex<double>* designSiteValues );

//Function sets up GPU for Gaussian Correlation calculations between two matrices
cudaError_t calculateGaussianCorrelation(std::complex<double>* outputMatrix, std::complex<double> variance, std::complex<double> a, std::complex<double> theta, std::complex<double>* inMatrixOne, std::complex<double>* inMatrixTwo);

//Function sets up GPU for calculating the distance between two matrices
cudaError_t calculateDistanceBetweenMatrices(std::complex<double>* outputMatrix, std::complex<double>* inMatrixOne, std::complex<double>* inMatrixTwo);

//Function sets up GPU for calculating the weights vector
cudaError_t calculateWeightVector(std::complex<double>* outputVector, std::complex<double>* invertedCovarianceMatrix, std::complex<double>*);

//Function sets up GPU for calculating the inverse of a matrix using Gauss-Jordan elimination
cudaError_t invertMatrix(std::complex<double>* outputMatrix, std::complex<double>* inputMatrix, int dimension);

//Function sets the given matrix to the identity matrix
void createIdentityMatrix(std::complex<double>* matrix);

//Begin GPU Function Definitions
//CUDA function to calculate distance between two matrices
__global__ void calcDistanceBetMat(std::complex<double>* outMat, std::complex<double>* inMatOne, std::complex<double>* inMatTwo, int dimension);

//CUDA function to calculate the Gaussian Correlation of two matrices
__global__ void calculateGaussCorr(std::complex<double>* outMat, std::complex<double>* inMatOne, std::complex<double>* inMatTwo);

//CUDA functions to calculate the inverse of a given matrix using Gauss-Jordan elimination method
__global__ void nonDiagNormalize(std::complex<double>* inputMatrix, std::complex<double>* identityMatrix, int n, int i);
__global__ void diagNormalize(std::complex<double>* inputMatrix, std::complex<double>* identityMatrix, int n, int i);
__global__ void gaussJordan(std::complex<double>* inputMatrix, std::complex<double>* identityMatrix, int n, int i);
__global__ void setZero(std::complex<double>* inputMatrix, std::complex<double>* identityMatrix, int n, int i);

//Begin Function Implementations
cudaError_t metamodelSetup(std::complex<double>* outputValue, int dimension, std::complex<double> theta, std::complex<double> variance, std::complex<double> a, std::complex<double>* designSite, std::complex<double>* testSite, std::complex<double>* designSiteValues) {
    
    //Begin variable definitions for data to be passed to GPU
    std::complex<double>* deviceOutVal = 0;
    std::complex<double>* deviceDesignSite = 0;
    std::complex<double>* deviceDesignSiteValues = 0;
    std::complex<double>* deviceTestSite = 0;
    std::complex<double>* deviceIdentityMatrix = 0;
    cudaError_t cudaStatus = cudaSuccess;
    int matrixMemoryAllocationSize = pow(dimension, 2);
    int vectorMemoryAllocationSize = dimension;

    //Create a dynamic allocation of memory for the identity matrix
    std::complex<double>* identityMatrix = (std::complex<double>*) malloc(matrixMemoryAllocationSize * sizeof(std::complex<double>));

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess){
        goto Error;
    }

    //Begin GPU memory allocation
    cudaStatus = cudaMalloc((void**)&deviceIdentityMatrix, matrixMemoryAllocationSize * sizeof(std::complex<double>));
    if (cudaStatus != cudaSuccess) {
        goto Error;
    }



    //Define error state
Error:
    cudaFree(deviceOutVal);
    cudaFree(deviceDesignSite);
    cudaFree(deviceDesignSiteValues);
    cudaFree(deviceIdentityMatrix);
    free(identityMatrix);

    return cudaStatus;
}

//Begin CUDA Function Implementations
__global__ void nonDiagNormalize(std::complex<double>* inputMatrix, std::complex<double>* identityMatrix, int n, int i) {
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

__global__ void diagNormalize(std::complex<double>* inputMatrix, std::complex<double>* identityMatrix, int n, int i) {
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

__global__ void gaussJordan(std::complex<double>* inputMatrix, std::complex<double>* identityMatrix, int n, int i) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if ((x < n) && (y < n)) {
        if (x != i) {
            identityMatrix[x * n + y] -= identityMatrix[i * n + y] * inputMatrix[x * n + i];
            if (y != i) {
                inputMatrix[x * n + y] -= inputMatrix[i * n + y] * inputMatrix[x * n + i];
            }
        }
    }

    return;
}

__global__ void setZero(std::complex<double>* inputMatrix, std::complex<double>* identityMatrix, int n, int i) {
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

//Old functions
/*
//Function calculates the average of the given matrix
double calculateAverage(double inputMatrix[], int dimension);

//Function prints the given matrix on the GPU
__global__ void printWithGPU(std::complex<double>* matrix, int dimension);

//Function calculates the covariance of the given matrices and places them into a new output matrix
__global__ void calculateCovariance(std::complex<double>* outputMatrix, std::complex<double>* X, double meanX, std::complex<double>* Y, double meanY, int dimension);

//Function calculates the Kriging Estimator for the Global Mean
__global__ void calculateKrigingEstimator(std::complex<double>* result, std::complex<double>* covarianceMatrix, std::complex<double>* sampleMatrix);

//Functions calculate the dot product of two vectors
__global__ void calculateDotProduct(std::complex<double>* result, std::complex<double>* inputOne, int xOne, int yOne, std::complex<double>* inputTwo, int xTwo, int yTwo);

// Helper function for using CUDA to add vectors in parallel.
cudaError_t cudaSetup(std::complex<double>* c, std::complex<double>* a, std::complex<double>* b, std::complex<double>* identityMatrix, unsigned int size)
{
    std::complex<double>* deviceFirstInput = 0;
    std::complex<double>* deviceSecondInput = 0;
    std::complex<double>* deviceOutput = 0;
    std::complex<double>* deviceIdentityMatrix = 0;
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

__global__ void printWithGPU(std::complex<double>* matrix, int dimension) {

    printf("Hello from <%d, %d>: %f\n", blockIdx.x, threadIdx.x, matrix[blockIdx.x + (threadIdx.x * dimension)]);
    return;
}

__global__ void calculateCovariance(std::complex<double>* outputMatrix, std::complex<double>* X, double meanX, std::complex<double>* Y, double meanY, int dimension) {
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

__global__ void calculateKrigingEstimator(std::complex<double>* result, std::complex<double>* covarianceMatrix, std::complex<double>* sampleMatrix) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;



    return;
}

*/