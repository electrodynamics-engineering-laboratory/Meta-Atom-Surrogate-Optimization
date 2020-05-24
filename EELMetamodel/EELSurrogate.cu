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
cudaError_t calculateGaussianCorrelation(std::complex<double>* outputMatrix, std::complex<double>* inMatrix, std::complex<double> variance, std::complex<double> a, std::complex<double> theta, int dimension);

//Function sets up GPU for calculating the distance between two matrices
cudaError_t calculateDistanceBetweenMatrices(std::complex<double>* outputMatrix, std::complex<double>* inMatrixOne, std::complex<double>* inMatrixTwo, int dimension);

//Function sets up GPU for calculating the distance between a matrix and a vector
cudaError_t calculateDistanceBetweenMatrixVector(std::complex<double>* outputMatrix, std::complex<double>* inMatrix, std::complex<double>* inVector, int dimension);

//Function sets up GPU for calculating the weights vector
cudaError_t calculateWeightVector(std::complex<double>* outputVector, std::complex<double>* invertedCovarianceMatrix, std::complex<double>* covarianceVector, int dimension);

//Function sets up GPU for extending a matrix
cudaError_t extendMatrix(std::complex<double>* outputMatrix, std::complex<double>* inputMatrix, int dimension);

//Function sets up GPU for calculating the inverse of a matrix using Gauss-Jordan elimination
cudaError_t invertMatrix(std::complex<double>* outputMatrix, std::complex<double>* inputMatrix, int dimension);

//Function sets the given matrix to the identity matrix
cudaError_t createIdentityMatrix(std::complex<double>* matrix, int dimension);

//Begin GPU Function Definitions
//CUDA function to calculate distance between two matrices
__global__ void calcDistanceBetMats(std::complex<double>* outMat, std::complex<double>* inMatOne, std::complex<double>* inMatTwo, int dimension);

//CUDA function to calculate the distance between a matrix and a vector
__global__ void calcDistanceBetMatVec(std::complex<double>* outMat, std::complex<double>* inMat, std::complex<double>* inVec, int dimension);

//CUDA function to calculate the Gaussian Correlation of two matrices
__global__ void calcGaussCorr(std::complex<double>* outMat, std::complex<double>* inMat, int dimension, std::complex<double> variance, std::complex<double> a, std::complex<double> theta);

//CUDA function to normalize a given matrix to a value
__global__ void normalizeMatrix(std::complex<double>* outMat, std::complex<double>* inMat, std::complex<double> normalizingValue, int dimension);

//CUDA function to extend a given matrix
__global__ void extendMat(std::complex<double>* outMat, std::complex<double>* inMat, int dimension);

//CUDA functions to calculate the inverse of a given matrix using Gauss-Jordan elimination method
__global__ void nonDiagNormalize(std::complex<double>* inputMatrix, std::complex<double>* identityMatrix, int n, int i);
__global__ void diagNormalize(std::complex<double>* inputMatrix, std::complex<double>* identityMatrix, int n, int i);
__global__ void gaussJordan(std::complex<double>* inputMatrix, std::complex<double>* identityMatrix, int n, int i);
__global__ void setZero(std::complex<double>* inputMatrix, std::complex<double>* identityMatrix, int n, int i);

//CUDA function to calculate the dot product of two inputs
__global__ void multiplyMatrix(std::complex<double>* output, std::complex<double>* firInput, std::complex<double>* secInput, int dimension);

//CUDA function to create an identity matrix of the given dimension
__global__ void createIdentMat(std::complex<double>* matrix, int dimension);

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

    //Create a dynamic allocation of memory for the identity matrix and populate with values
    std::complex<double>* identityMatrix = (std::complex<double>*) malloc(matrixMemoryAllocationSize * sizeof(std::complex<double>));
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            if (i == j) {
                identityMatrix[i + j * dimension] = 1;
            }
            else {
                identityMatrix[i + j * dimension] = 0;
            }
        }
    }

    //Create a dynamic allocation of memory for a temporary holding matrix
    std::complex<double>* tempMatrix = (std::complex<double>*) malloc(matrixMemoryAllocationSize * sizeof(std::complex<double>));
    for (int i = 0; i < pow(dimension,2); i++) {
        tempMatrix[i] = 0;
    }

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess){
        goto SetupError;
    }

    //Begin GPU memory allocation

    //Allocate GPU memory for the identity matrix
    cudaStatus = cudaMalloc((void**)&deviceIdentityMatrix, matrixMemoryAllocationSize * sizeof(std::complex<double>));
    if (cudaStatus != cudaSuccess) {
        goto SetupError;
    }

    //Allocate GPU memory for the test site matrix
    cudaStatus = cudaMalloc((void**)&deviceTestSite, matrixMemoryAllocationSize * sizeof(std::complex<double>));
    if (cudaStatus != cudaSuccess) {
        goto SetupError;
    }

    //Allocate GPU memory for the design site values
    cudaStatus = cudaMalloc((void**)&deviceDesignSiteValues, matrixMemoryAllocationSize * sizeof(std::complex<double>));
    if (cudaStatus != cudaSuccess) {
        goto SetupError;
    }

    //Allocate GPU memory for the design site matrix
    cudaStatus = cudaMalloc((void**)&deviceDesignSite, matrixMemoryAllocationSize * sizeof(std::complex<double>));
    if (cudaStatus != cudaSuccess) {
        goto SetupError;
    }

    //Allocate GPU memory for the output matrix
    cudaStatus = cudaMalloc((void**)&deviceOutVal, matrixMemoryAllocationSize * sizeof(std::complex<double>));
    if (cudaStatus != cudaSuccess) {
        goto SetupError;
    }

    //Define error state
SetupError:
    cudaFree(deviceOutVal);
    cudaFree(deviceDesignSite);
    cudaFree(deviceDesignSiteValues);
    cudaFree(deviceIdentityMatrix);
    free(identityMatrix);

    return cudaStatus;
}

cudaError_t calculateGaussianCorrelation(std::complex<double>* outputMatrix, std::complex<double>* inMatrix, std::complex<double> variance, std::complex<double> a, std::complex<double> theta, int dimension) {
    std::complex<double>* deviceOutVal = 0;
    std::complex<double>* deviceInMat = 0;
    cudaError_t cudaStatus = cudaSuccess;
    int matrixMemoryAllocationSize = pow(dimension, 2);
    int vectorMemoryAllocationSize = dimension;

    cudaStatus = cudaMalloc((void**)&deviceOutVal, matrixMemoryAllocationSize * sizeof(std::complex<double>));
    if (cudaStatus != cudaSuccess) {
        goto CorrError;
    }

    cudaStatus = cudaMalloc((void**)&deviceInMat, matrixMemoryAllocationSize * sizeof(std::complex<double>));
    if (cudaStatus != cudaSuccess) {
        goto CorrError;
    }

    cudaStatus = cudaMemcpy(deviceInMat, inMatrix, matrixMemoryAllocationSize * sizeof(std::complex<double>), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        goto CorrError;
    }

    //Perform calculation on the GPU and catch any error
    calcGaussCorr << < dimension, dimension >> > (deviceOutMat, deviceInMatOne, deviceInMatTwo, dimension);
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        goto CorrError;
    }

    //Synchronize device to enable a copy of the result
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        goto CorrError;
    }

    cudaStatus = cudaMemcpy(outputMatrix, deviceOutVal, matrixMemoryAllocationSize * sizeof(std::complex<double>), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        goto CorrError;
    }

CorrError:
    cudaFree(deviceOutVal);
    cudaFree(deviceInMat);

    return cudaStatus;
}

cudaError_t calculateDistanceBetweenMatrices(std::complex<double>* outputMatrix, std::complex<double>* inMatrixOne, std::complex<double>* inMatrixTwo, int dimension) {
    std::complex<double>* deviceOutMat = 0;
    std::complex<double>* deviceInMatOne = 0;
    std::complex<double>* deviceInMatTwo = 0;
    cudaError_t cudaStatus = cudaSuccess;
    int matrixMemoryAllocationSize = pow(dimension, 2);

    //Begin allocation of memory on device
    cudaStatus = cudaMalloc((void**)&deviceOutMat, matrixMemoryAllocationSize * sizeof(std::complex<double>));
    if (cudaStatus != cudaSuccess) {
        goto MatDistError;
    }

    cudaStatus = cudaMalloc((void**)&deviceInMatOne, matrixMemoryAllocationSize * sizeof(std::complex<double>));
    if (cudaStatus != cudaSuccess) {
        goto MatDistError;
    }

    cudaStatus = cudaMalloc((void**)&deviceInMatTwo, matrixMemoryAllocationSize * sizeof(std::complex<double>));
    if (cudaStatus != cudaSuccess) {
        goto MatDistError;
    }

    //Begin copy of CPU local data to GPU
    cudaStatus = cudaMemcpy(deviceInMatOne, inMatrixOne, matrixMemoryAllocationSize * sizeof(std::complex<double>), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        goto MatDistError;
    }

    cudaStatus = cudaMemcpy(deviceInMatTwo, inMatrixTwo, matrixMemoryAllocationSize * sizeof(std::complex<double>), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        goto MatDistError;
    }

    //Perform calculation on the GPU and catch any error
    calcDistanceBetMats <<< dimension, dimension >>> (deviceOutMat, deviceInMatOne, deviceInMatTwo, dimension, twoMatrices);
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        goto MatDistError;
    }

    //Synchronize device to enable a copy of the result
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        goto MatDistError;
    }

    //Copy data from GPU address to CPU address
    cudaStatus = cudaMemcpy(outputMatrix, deviceOutMat, matrixMemoryAllocationSize * sizeof(std::complex<double>), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        goto MatDistError;
    }

MatDistError:
    cudaFree(deviceOutMat);
    cudaFree(deviceInMatOne);
    cudaFree(deviceInMatTwo);

    return cudaStatus;
}

cudaError_t calculateDistanceBetweenMatrixVector(std::complex<double>* outputMatrix, std::complex<double>* inMatrix, std::complex<double>* inVector, int dimension) {
    std::complex<double>* deviceOutMat = 0;
    std::complex<double>* deviceInMat = 0;
    std::complex<double>* deviceInVec = 0;
    cudaError_t cudaStatus = cudaSuccess;
    int matrixMemoryAllocationSize = pow(dimension, 2);

    //Begin allocation of memory on device
    cudaStatus = cudaMalloc((void**)&deviceOutMat, matrixMemoryAllocationSize * sizeof(std::complex<double>));
    if (cudaStatus != cudaSuccess) {
        goto VecDistError;
    }

    cudaStatus = cudaMalloc((void**)&deviceInMat, matrixMemoryAllocationSize * sizeof(std::complex<double>));
    if (cudaStatus != cudaSuccess) {
        goto VecDistError;
    }

    cudaStatus = cudaMalloc((void**)&deviceInVec, matrixMemoryAllocationSize * sizeof(std::complex<double>));
    if (cudaStatus != cudaSuccess) {
        goto VecDistError;
    }

    //Begin copy of CPU local data to GPU
    cudaStatus = cudaMemcpy(deviceInMat, inMatrix, matrixMemoryAllocationSize * sizeof(std::complex<double>), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        goto VecDistError;
    }

    cudaStatus = cudaMemcpy(deviceInVec, inVector, matrixMemoryAllocationSize * sizeof(std::complex<double>), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        goto VecDistError;
    }

    //Perform calculation on the GPU and catch any error
    calcDistanceBetMatVec <<< dimension, 1 >>> (deviceOutMat, deviceInMat, deviceInVec, dimension);
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        goto VecDistError;
    }

    //Synchronize device to enable a copy of the result
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        goto VecDistError;
    }

    //Copy data from GPU address to CPU address
    cudaStatus = cudaMemcpy(outputMatrix, deviceOutMat, matrixMemoryAllocationSize * sizeof(std::complex<double>), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        goto VecDistError;
    }

VecDistError:
    cudaFree(deviceOutMat);
    cudaFree(deviceInMat);
    cudaFree(deviceInVec);

    return cudaStatus;
}

cudaError_t calculateWeightVector(std::complex<double>* outputVectorMatrix, std::complex<double>* invertedCovarianceMatrix, std::complex<double>* covarianceVectorMatrix, int dimension) {
    std::complex<double>* deviceOutMat = 0;
    std::complex<double>* deviceInvCovMat = 0;
    std::complex<double>* deviceCovVecMat = 0;
    cudaError_t cudaStatus = cudaSuccess;
    int matrixMemoryAllocationSize = pow(dimension, 2);
    int vectorMemoryAllocationSize = dimension;

    //Begin GPU memory allocation
    cudaStatus = cudaMalloc((void**)&deviceOutMat, matrixMemoryAllocationSize * sizeof(std::complex<double>));
    if (cudaStatus != cudaSuccess) {
        goto WeightError;
    }

    cudaStatus = cudaMalloc((void**)&deviceInvCovMat, matrixMemoryAllocationSize * sizeof(std::complex<double>));
    if (cudaStatus != cudaSuccess) {
        goto WeightError;
    }

    cudaStatus = cudaMalloc((void**)&deviceCovVecMat, vectorMemoryAllocationSize * sizeof(std::complex<double>));
    if (cudaStatus != cudaSuccess) {
        goto WeightError;
    }

    //Copy data from local CPU to GPU
    cudaStatus = cudaMemcpy(deviceInvCovMat, invertedCovarianceMatrix, matrixMemoryAllocationSize*sizeof(std::complex<double>), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        goto WeightError;
    }

    cudaStatus = cudaMemcpy(deviceCovVecMat, covarianceVectorMatrix, vectorMemoryAllocationSize * sizeof(std::complex<double>), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        goto WeightError;
    }

    multiplyMatrix<<<dimension, dimension>>>(deviceOutMat, deviceInvCovMat, deviceCovVecMat)
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        goto WeightError;
    }

    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        goto WeightError;
    }

    cudaStatus = cudaMemcpy(outputVectorMatrix, deviceOutMat, matrixMemoryAllocationSize * sizeof(std::complex<double>), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        goto WeightError;
    }


WeightError:
    cudaFree(deviceOutMat);
    cudaFree(deviceInvCovMat);
    cudaFree(deviceCovVecMat);

    return cudaStatus;
}

cudaError_t extendMatrix(std::complex<double>* outputMatrix, std::complex<double>* inputMatrix, int dimension) {
    std::complex<double>* deviceOutMat = 0;
    std::complex<double>* deviceInMat = 0;
    cudaError_t cudaStatus = cudaSuccess;
    int matrixMemoryAllocationSize = pow(dimension, 2);

    //Begin allocation of memory on GPU device
    cudaStatus = cudaMalloc((void**)&deviceOutMat, matrixMemoryAllocationSize * sizeof(std::complex<double>));
    if (cudaStatus != cudaSuccess) {
        goto ExtendError;
    }

    cudaStatus = cudaMalloc((void**)&deviceInMat, matrixMemoryAllocationSize * sizeof(std::complex<double>));
    if (cudaStatus != cudaSuccess) {
        goto ExtendError;
    }

    //Begin copy of CPU local data to GPU
    cudaStatus = cudaMemcpy(deviceInMat, inputMatrix, matrixMemoryAllocationSize * sizeof(std::complex<double>), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        goto ExtendError;
    }

    //Perform calculation on the GPU and catch any error
    extendMat <<<dimension + 1, dimension + 1 >>> (deviceOutMat, deviceInMat, dimension);
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        goto ExtendError;
    }

    //Synchronize device to enable a copy of the result
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        goto ExtendError;
    }

    //Copy data from GPU address to CPU address
    cudaStatus = cudaMemcpy(outputMatrix, deviceOutMat, matrixMemoryAllocationSize * sizeof(std::complex<double>), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        goto ExtendError;
    }

ExtendError:
    cudaFree(deviceOutMat);
    cudaFree(deviceInMat);

    return cudaStatus;
}

cudaError_t invertMatrix(std::complex<double>* outputMatrix, std::complex<double>* inputMatrix, int dimension) {
    std::complex<double>* deviceInMat = 0;
    std::complex<double>* deviceIdenMat = 0;
    cudaError_t cudaStatus = cudaSuccess;
    int matrixMemoryAllocationSize = pow(dimension, 2);

    cudaStatus = cudaMalloc((void**)&deviceInMat, matrixMemoryAllocationSize * sizeof(std::complex<double>));
    if (cudaStatus != cudaSuccess) {
        goto InvertError;
    }

    cudaStatus = cudaMalloc((void**)&deviceIdenMat, matrixMemoryAllocationSize * sizeof(std::complex<double>));
    if (cudaStatus != cudaSuccess) {
        goto InvertError;
    }

    cudaStatus = cudaMemcpy(deviceInMat, inputMatrix, matrixMemoryAllocationSize * sizeof(std::complex<double>), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        goto InvertError;
    }

    cudaStatus = cudaMemcpy(deviceIdenMat, outputMatrix, matrixMemoryAllocationSize * sizeof(std::complex<double>), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        goto InvertError;
    }

    for (int i = 0; i < dimension; i++) {
        nonDiagNormalize<<<dimension, dimension>>> (deviceInMat, deviceIdenMat, dimension, i);
        diagNormalize <<<dimension, dimension >>> (deviceInMat, deviceIdenMat, dimension, i);
        gaussJordan <<<dimension, dimension >>> (deviceInMat, deviceIdenMat, dimension, i);
        setZero <<<dimension, dimension >>> (deviceInMat, deviceIdenMat, dimension, i);
    }
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        goto InvertError;
    }

    cudaStatus = cudaMemcpy(outputMatrix, deviceIdenMat, matrixMemoryAllocationSize * sizeof(std::complex<double>), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        goto InvertError;
    }

    cudaStatus = cudaMemcpy(inputMatrix, deviceInMat, matrixMemoryAllocationSize * sizeof(std::complex<double>), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        goto InvertError;
    }

InvertError:
    cudaFree(deviceIdenMat);
    cudaFree(deviceInMat);

    return cudaStatus;
}

cudaError_t createIdentityMatrix(std::complex<double>* matrix, int dimension) {
    std::complex<double>* deviceMat = 0;
    cudaError_t cudaStatus = cudaSuccess;
    int matrixMemoryAllocationSize = pow(dimension, 2);

    cudaStatus = cudaMalloc((void**)&deviceMat, matrixMemoryAllocationSize * sizeof(std::complex<double>));
    if (cudaStatus != cudaSuccess) {
        goto IdentityError;
    }

    cudaStatus = cudaMemcpy(deviceMat, matrix, matrixMemoryAllocationSize * sizeof(std::complex<double>), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        goto IdentityError;
    }

    createIdentMat <<<dimension, dimension >> > (deviceMat, dimension);
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        goto IdentityError;
    }

    cudaStatus = cudaMemcpy(matrix, deviceMat, matrixMemoryAllocationSize * sizeof(std::complex<double>), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
    }


IdentityError:
    cudaFree(deviceMat);

    return cudaStatus;
}

//Begin CUDA Function Implementations
__global__ void calcDistanceBetMats(std::complex<double>* outMat, std::complex<double>* inMatOne, std::complex<double>* inMatTwo, int dimension) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    outMat[i + j * dimension] = std::pow(std::abs(inMatOne[i + j * dimension] - inMatTwo[i + j * dimension]),2);
    return;
}

__global__ void calcDistanceBetMatVec(std::complex<double>* outMat, std::complex<double>* inMat, std::complex<double>* inVec, int dimension) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = 0;
    
    outMat[i + j * dimension] = std::pow(std::abs(inMat[i + j * dimension] - inVec[i + j * dimension]), 2);
    return;
}

__global__ void calcGaussCorr(std::complex<double>* outMat, std::complex<double>* inMat, int dimension, std::complex<double> variance, std::complex<double> a, std::complex<double> theta) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    std::complex<double> negOne(-1, 0);

    outMat[i + j * dimension] = (variance - a) * std::exp(negOne * theta * inMat[i + j * dimension]);

    return;
}

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

__global__ void multiplyMatrix(std::complex<double>* output, std::complex<double>* firInput, std::complex<double>* secInput, int dimension) {
    int x = blockIdx.x;
    int y = threadIdx.x;

    //Multiple each element pair and then sum them together
    for (int i = 0; i < dimension; i++) {
        output[x + y * dimension] += firInput[i + y * dimension] * secInput[x + i * dimension];
    }

    return;
}

__global__ void extendMat(std::complex<double>* outMat, std::complex<double>* inMat, int dimension) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    //If both i & j equal dimension, the current index is the corner which should be set to zero
    if (i == dimension && j == dimension) {
        outMat[i + j * (dimension + 1)] = 0;
    }
    //If i or j are equal to dimension, the current index is the right or bottom edge of the matrix which should be set to one
    else if (i == dimension || j == dimension) {
        outMat[i + j * (dimension + 1)] = 1;
    }
    //Otherwise, simply copy the value over to the output matrix
    else {
        outMat[i + j * (dimension + 1)] = inMat[i + j * (dimension)];
    }
    return;
}

__global__ void normalizeMatrix(std::complex<double>* outMat, std::complex<double>* inMat, std::complex<double> normalizingValue, int dimension) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    outMat[i + j * dimension] = inMat[i + j * dimension] / normalizingValue;

    return;
}

__global__ void createIdentMat(std::complex<double>* matrix, int dimension) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i == j) {
        matrix[i + j * dimension] = 1;
    }
    else {
        matrix[i + j * dimension] = 0;
    }

    return;
}
