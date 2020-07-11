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
#include <math.h>
#include <cuda_fp16.h>
//#include <cuComplex.h> //Unable to find proper documentation for cuComplex functionality. As such, some of the documentation relating to complex values is meaningless. 

//Begin CPU Function Definitions
/*
Name: metamodelSetup
Description: Function sets up and performs calculations on the GPU that relate to a metamodel
Input:  outputValue - (CURRENTLY UNKNOWN AND WILL CHANGE) a pointer to a value that will be used for output 
        dimension - an integer value representing the number of rows for all matrices and vectors and the number columns for all matrices
        theta - (THIS MIGHT BE AN INCORRECT DESCRIPTION) a value in radians that represents the angle between the two matrices, this value should only ever be real
        variance - the value of the variance between 
        designSite - a dimension-by-dimension matrix that represents the datapoints previously collected to be used to estimate a new value based on test site
        testSite - a dimension-by-dimension matrix with only the first column that make up the vector that represents the points at which a new value is to be calculated
        designSiteValues - a dimension-by-dimension matrix that represent the values at particular design sites. Each matrix has a 1-to-1 correlation of site to value (i.e. designSite[0] is paired with designSiteValues[0])
Output: a value that is any error that the GPU experienced upon attempting calculations
*/
//Function sets up and performs calculations on the GPU
double metamodelSetup(int dimension, double theta, double variance, double a, double* designSite, double* testSite, double* designSiteValues );

//Function sets up GPU for covariance calculations between two matrices through Gaussian Correlation 
cudaError_t calculateGaussianCorrelation(double* outputMatrix, double* inMatrix, double variance, double a, double theta, int dimension);

//Function sets up GPU for calculating the distance between two matrices
cudaError_t calculateDistanceBetweenMatrices(double* outputMatrix, double* inMatrixOne, double* inMatrixTwo, int dimension);

//Function sets up GPU for calculating the distance between a matrix and a vector
cudaError_t calculateDistanceBetweenMatrixVector(double* outputMatrix, double* inMatrix, double* inVector, int dimension);

//Function sets up GPU for calculating the weights vector
cudaError_t calculateWeightVector(double* outputVector, double* invertedCovarianceMatrix, double* covarianceVector, int dimension);

//Function sets up GPU for extending a matrix
cudaError_t extendMatrix(double* outputMatrix, double* inputMatrix, int dimension);

//Function sets up GPU for calculating the inverse of a matrix using Gauss-Jordan elimination
cudaError_t invertMatrix(double* outputMatrix, double* inputMatrix, int dimension);

//Function sets the given matrix to the identity matrix
cudaError_t createIdentityMatrix(double* matrix, int dimension);

//Function multiples two matrices
cudaError_t multiplyMatrices(double* output, double* inputMatrix, int dimension);

//Begin GPU Function Definitions
//CUDA function to calculate distance between t6wo matrices
__global__ void calcDistanceBetMats(double* outMat, double* inMatOne, double* inMatTwo, int dimension);

//CUDA function to calculate the distance between a matrix and a vector
__global__ void calcDistanceBetMatVec(double* outMat, double* inMat, double* inVec, int dimension);

//CUDA function to calculate the Gaussian Correlation of two matrices
__global__ void calcGaussCorr(double* outMat, double* inMat, int dimension, double variance, double a, double theta);

//CUDA function to normalize a given matrix to a value
__global__ void normalizeMatrix(double* outMat, double* inMat, double normalizingValue, int dimension);

//CUDA function to extend a given matrix
__global__ void extendMat(double* outMat, double* inMat, int dimension);

//CUDA functions to calculate the inverse of a given matrix using Gauss-Jordan elimination method
__global__ void nonDiagNormalize(double* inputMatrix, double* identityMatrix, int n, int i);
__global__ void diagNormalize(double* inputMatrix, double* identityMatrix, int n, int i);
__global__ void gaussJordan(double* inputMatrix, double* identityMatrix, int n, int i);
__global__ void setZero(double* inputMatrix, double* identityMatrix, int n, int i);

//CUDA function to calculate the dot product of two inputs
__global__ void multiplyMatrix(double* output, double* firInput, double* secInput, int dimension);

//CUDA function to create an identity matrix of the given dimension
__global__ void createIdentMat(double* matrix, int dimension);

//Testing function to print out arrays after computation.
void printMatrix(double inArray[], int dimension);

//Begin Function Implementations
double metamodelSetup(int dimension, double theta, double variance, double a, double* designSite, double* testSite, double* designSiteValues) {
    
    //Begin variable definitions for data to be passed to GPU
    cudaError_t cudaStatus = cudaSuccess;
    int matrixMemoryAllocationSize = pow(dimension+1, 2);
    int extendedMatrixMemoryAllocationSize = pow(dimension+1, 2);
    int vectorMemoryAllocationSize = dimension+1;
    int extendedVectorMemoryAllocationSize = dimension + 1;
    double outputValue = 0;

    //Create a dynamic allocation of memory for the identity matrix and populate with values
    double* identityMatrix = (double*) malloc(matrixMemoryAllocationSize * sizeof(double));
    
    //Create a dynamic allocation of memory for a temporary holding matrix
    double* tempMatrixOne = (double*) malloc(extendedMatrixMemoryAllocationSize * sizeof(double));
    double* tempMatrixTwo = (double*) malloc(extendedMatrixMemoryAllocationSize * sizeof(double));
    for (int i = 0; i < pow(dimension+1,2); i++) {
        tempMatrixOne[i] = 0;
        tempMatrixTwo[i] = 0;
    }

    //std::cout << "SETUP: Blank matrices." << std::endl;
    //printMatrix(tempMatrixOne, dimension);
    //printMatrix(tempMatrixTwo, dimension);

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess){
        goto SetupError;
    }

    //Create identity matrix on GPU, input is also the output
    cudaStatus = createIdentityMatrix(identityMatrix, dimension);
    if (cudaStatus != cudaSuccess) {
        goto SetupError;
    }

    //std::cout << "SETUP: Identity matrix." << std::endl;
    //printMatrix(identityMatrix, dimension);

    //Calculate distance between design sites and values at design sites. tempMatrixOne will hold the output
    cudaStatus = calculateDistanceBetweenMatrices(tempMatrixOne, designSite, designSiteValues, dimension);
    if (cudaStatus != cudaSuccess) {
        goto SetupError;
    }

    //printf("SETUP: Distance between designSite and designSiteValues\n");
    //printMatrix(tempMatrixOne, dimension);

    //Calculate distance between test site and design site. tempMatrixTwo will hold the output
    cudaStatus = calculateDistanceBetweenMatrixVector(tempMatrixTwo, designSite, testSite, dimension);
    if (cudaStatus != cudaSuccess) {
        goto SetupError;
    }
    
    //printf("SETUP: Distance between designSite and testSite\n");
    //printMatrix(tempMatrixTwo, dimension);

   // printf("SETUP: Dimension %d\n", dimension);
    //Calculate the covariance between design sites and values at design sites. tempMatrixOne will hold the output
    cudaStatus = calculateGaussianCorrelation(tempMatrixOne, tempMatrixOne, variance, a, theta, dimension);
    if (cudaStatus != cudaSuccess) {
        goto SetupError;
    }

    printf("SETUP: Matrix Gaussian Correlation\n");
    printMatrix(tempMatrixOne, dimension);

    //Calculate the covariance between test sites and design sites. tempMatrixTwo will hold the output
    cudaStatus = calculateGaussianCorrelation(tempMatrixTwo, tempMatrixTwo, variance, a, theta, dimension);
    if (cudaStatus != cudaSuccess) {
        goto SetupError;
    }

    //printf("SETUP: Vector Gaussian Correlation\n");
    //printMatrix(tempMatrixTwo, dimension);

    //Extend the covariance matrix between the design site and design site values
    cudaStatus = extendMatrix(tempMatrixOne, tempMatrixOne, dimension);
    if (cudaStatus != cudaSuccess) {
        goto SetupError;
    }

    printf("SETUP: Extend Matrix\n");
    printMatrix(tempMatrixOne, dimension+1);

    //Extend the covariance vector between test site and design sites. 
    tempMatrixTwo[dimension + 0*dimension] = 1; //Add 1 to the last row of the matrix, unclear if this is correct
    
    printf("SETUP: Extend Vector\n");
    printMatrix(tempMatrixTwo, dimension+1);

    //Calculate inverse of extended covariance matrix
    //cudaStatus = invertMatrix(tempMatrixOne, tempMatrixOne, dimension + 1);
    if (cudaStatus != cudaSuccess) {
        goto SetupError;
    }

    //Calculate extended weights vector, tempMatrixTwo will hold the result
    //cudaStatus = calculateWeightVector(tempMatrixTwo, tempMatrixOne, tempMatrixTwo, dimension + 1);
    if (cudaStatus != cudaSuccess) {
        goto SetupError;
    }

    //Calculate estimate value at test site, the ultimate output. 
    //cudaStatus = multiplyMatrices(tempMatrixTwo, designSiteValues, dimension); //Only consider elements within the dimension as the final value in the weight matrix is the lamda value (AKA not needed)
    if (cudaStatus != cudaSuccess) {
        goto SetupError;
    }

    //Grab the first, and only, value of tempMatrixTwo as the ultimate output value
    //outputValue = tempMatrixTwo[0];

    //Define error state
    //Need to manage/report error state so that the output value can be returned NOT an error status
SetupError:
    if (cudaStatus != cudaSuccess) {
        std::cout << "SETUP: Device failed" << std::endl;
        std::cout << "SETUP: CUDA Error Code -> " << cudaGetErrorString(cudaStatus) << std::endl;
    }
    free(identityMatrix);
    free(tempMatrixOne);
    free(tempMatrixTwo);
    
    return outputValue;
}

//Use the outputMatrix is initially an input, but the matrix is then overwritten for the output to save memory. Unclear if this is entirely necessary to prevent overflowing memory.
cudaError_t calculateGaussianCorrelation(double* outputMatrix, double* inMatrix, double variance, double a, double theta, int dimension) {
    
    //Allocate two pointers that will be used by the GPU for calculations.
    double* deviceOutMat = 0;
    double* deviceInMat = 0;

    //Create status variable for errors
    cudaError_t cudaStatus = cudaSuccess;

    //Create two integers that represent size to allocate on the GPU for matrices and vectors
    int matrixMemoryAllocationSize = pow(dimension, 2);
    int vectorMemoryAllocationSize = dimension;

    //Allocate memory on the GPU for the appropriate matrices
    cudaStatus = cudaMalloc((void**)&deviceOutMat, matrixMemoryAllocationSize * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        goto CorrError;
    }

    cudaStatus = cudaMalloc((void**)&deviceInMat, matrixMemoryAllocationSize * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        goto CorrError;
    }
    
    //Copy the CPU matrices to the GPU to allow for calculations
    cudaStatus = cudaMemcpy(deviceInMat, inMatrix, matrixMemoryAllocationSize * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        goto CorrError;
    }

    cudaStatus = cudaMemcpy(deviceOutMat, outputMatrix, matrixMemoryAllocationSize * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        goto CorrError;
    }
    
    //Call the GPU function with the appropriate number of blocks and threads to perform calculation on the GPU. Catch any error returned.
    calcGaussCorr <<< dimension, dimension >>> (deviceOutMat, deviceInMat, dimension, variance, a, theta);
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        goto CorrError;
    }

    //Synchronize device to enable a copy of the result
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        goto CorrError;
    }

    //Copy data from the GPU to the CPU matrix.
    cudaStatus = cudaMemcpy(outputMatrix, deviceOutMat, matrixMemoryAllocationSize * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        goto CorrError;
    }

CorrError:
    cudaFree(deviceOutMat);
    cudaFree(deviceInMat);

    return cudaStatus;
}

cudaError_t calculateDistanceBetweenMatrices(double* outputMatrix, double* inMatrixOne, double* inMatrixTwo, int dimension) {
    double* deviceOutMat = 0;
    double* deviceInMatOne = 0;
    double* deviceInMatTwo = 0;
    int* deviceDimension = 0;
    cudaError_t cudaStatus = cudaSuccess;
    int matrixMemoryAllocationSize = pow(dimension, 2);

    //Begin allocation of memory on device
    cudaStatus = cudaMalloc((void**)&deviceOutMat, matrixMemoryAllocationSize * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        goto MatDistError;
    }

    cudaStatus = cudaMalloc((void**)&deviceInMatOne, matrixMemoryAllocationSize * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        goto MatDistError;
    }

    cudaStatus = cudaMalloc((void**)&deviceInMatTwo, matrixMemoryAllocationSize * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        goto MatDistError;
    }

    //Begin copy of CPU local data to GPU
    cudaStatus = cudaMemcpy(deviceInMatOne, inMatrixOne, matrixMemoryAllocationSize * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        goto MatDistError;
    }

    cudaStatus = cudaMemcpy(deviceInMatTwo, inMatrixTwo, matrixMemoryAllocationSize * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        goto MatDistError;
    }

    //Perform calculation on the GPU and catch any error
    calcDistanceBetMats <<< dimension, dimension >>> (deviceOutMat, deviceInMatOne, deviceInMatTwo, dimension);
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
    cudaStatus = cudaMemcpy(outputMatrix, deviceOutMat, matrixMemoryAllocationSize * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        goto MatDistError;
    }

MatDistError:
    cudaFree(deviceOutMat);
    cudaFree(deviceInMatOne);
    cudaFree(deviceInMatTwo);

    return cudaStatus;
}

cudaError_t calculateDistanceBetweenMatrixVector(double* outputMatrix, double* inMatrix, double* inVector, int dimension) {
    double* deviceOutMat = 0;
    double* deviceInMat = 0;
    double* deviceInVec = 0;
    cudaError_t cudaStatus = cudaSuccess;
    int matrixMemoryAllocationSize = pow(dimension, 2);

    //Begin allocation of memory on device
    cudaStatus = cudaMalloc((void**)&deviceOutMat, matrixMemoryAllocationSize * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        goto VecDistError;
    }

    cudaStatus = cudaMalloc((void**)&deviceInMat, matrixMemoryAllocationSize * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        goto VecDistError;
    }

    cudaStatus = cudaMalloc((void**)&deviceInVec, matrixMemoryAllocationSize * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        goto VecDistError;
    }

    //Begin copy of CPU local data to GPU
    cudaStatus = cudaMemcpy(deviceInMat, inMatrix, matrixMemoryAllocationSize * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        goto VecDistError;
    }

    cudaStatus = cudaMemcpy(deviceInVec, inVector, matrixMemoryAllocationSize * sizeof(double), cudaMemcpyHostToDevice);
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
    cudaStatus = cudaMemcpy(outputMatrix, deviceOutMat, matrixMemoryAllocationSize * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        goto VecDistError;
    }

VecDistError:
    cudaFree(deviceOutMat);
    cudaFree(deviceInMat);
    cudaFree(deviceInVec);

    return cudaStatus;
}

cudaError_t calculateWeightVector(double* outputVectorMatrix, double* invertedCovarianceMatrix, double* covarianceVectorMatrix, int dimension) {
    double* deviceOutMat = 0;
    double* deviceInvCovMat = 0;
    double* deviceCovVecMat = 0;
    cudaError_t cudaStatus = cudaSuccess;
    int matrixMemoryAllocationSize = pow(dimension, 2);
    int vectorMemoryAllocationSize = dimension;

    //Begin GPU memory allocation
    cudaStatus = cudaMalloc((void**)&deviceOutMat, matrixMemoryAllocationSize * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        goto WeightError;
    }

    cudaStatus = cudaMalloc((void**)&deviceInvCovMat, matrixMemoryAllocationSize * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        goto WeightError;
    }

    cudaStatus = cudaMalloc((void**)&deviceCovVecMat, vectorMemoryAllocationSize * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        goto WeightError;
    }

    //Copy data from local CPU to GPU
    cudaStatus = cudaMemcpy(deviceInvCovMat, invertedCovarianceMatrix, matrixMemoryAllocationSize*sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        goto WeightError;
    }

    cudaStatus = cudaMemcpy(deviceCovVecMat, covarianceVectorMatrix, vectorMemoryAllocationSize * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        goto WeightError;
    }

    multiplyMatrix << <dimension, dimension >> > (deviceOutMat, deviceInvCovMat, deviceCovVecMat, dimension);
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        goto WeightError;
    }

    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        goto WeightError;
    }

    cudaStatus = cudaMemcpy(outputVectorMatrix, deviceOutMat, matrixMemoryAllocationSize * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        goto WeightError;
    }


WeightError:
    cudaFree(deviceOutMat);
    cudaFree(deviceInvCovMat);
    cudaFree(deviceCovVecMat);

    return cudaStatus;
}

cudaError_t extendMatrix(double* outputMatrix, double* inputMatrix, int dimension) {
    double* deviceOutMat = 0;
    double* deviceInMat = 0;
    cudaError_t cudaStatus = cudaSuccess;
    int matrixMemoryAllocationSize = pow(dimension+1, 2);

    //Swap matrix locations to prepare for "extension" of the matrix. 
    for (int i = dimension; i >= 0; i--) {
        for (int j = dimension; j >= 0; j--) {
            //printf("OLD [%d], NEW [%d]\n", (i + j * dimension), (i + j * (dimension + 1)));
            inputMatrix[i + j * (dimension + 1)] = inputMatrix[i + j * dimension];
        }
    }

    //Begin allocation of memory on GPU device
    cudaStatus = cudaMalloc((void**)&deviceOutMat, matrixMemoryAllocationSize * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        goto ExtendError;
    }

    cudaStatus = cudaMalloc((void**)&deviceInMat, matrixMemoryAllocationSize * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        goto ExtendError;
    }

    //Begin copy of CPU local data to GPU
    cudaStatus = cudaMemcpy(deviceInMat, inputMatrix, matrixMemoryAllocationSize * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        goto ExtendError;
    }

    //Perform calculation on the GPU and catch any error
    extendMat <<<dimension+ 1 , dimension+ 1>>> (deviceOutMat, deviceInMat, dimension+1); 
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
    printf("EXTEND: BEFORE\n");
    printMatrix(outputMatrix, dimension + 1);
    cudaStatus = cudaMemcpy(outputMatrix, deviceOutMat, matrixMemoryAllocationSize * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        goto ExtendError;
    }

    printf("EXTEND: AFTER\n");
    printMatrix(outputMatrix, dimension + 1);
ExtendError:
    cudaFree(deviceOutMat);
    cudaFree(deviceInMat);

    return cudaStatus;
}

cudaError_t invertMatrix(double* outputMatrix, double* inputMatrix, int dimension) {
    double* deviceInMat = 0;
    double* deviceIdenMat = 0;
    cudaError_t cudaStatus = cudaSuccess;
    int matrixMemoryAllocationSize = pow(dimension, 2);

    cudaStatus = cudaMalloc((void**)&deviceInMat, matrixMemoryAllocationSize * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        goto InvertError;
    }

    cudaStatus = cudaMalloc((void**)&deviceIdenMat, matrixMemoryAllocationSize * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        goto InvertError;
    }

    cudaStatus = cudaMemcpy(deviceInMat, inputMatrix, matrixMemoryAllocationSize * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        goto InvertError;
    }

    cudaStatus = cudaMemcpy(deviceIdenMat, outputMatrix, matrixMemoryAllocationSize * sizeof(double), cudaMemcpyHostToDevice);
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

    cudaStatus = cudaMemcpy(outputMatrix, deviceIdenMat, matrixMemoryAllocationSize * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        goto InvertError;
    }

    cudaStatus = cudaMemcpy(inputMatrix, deviceInMat, matrixMemoryAllocationSize * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        goto InvertError;
    }

InvertError:
    cudaFree(deviceIdenMat);
    cudaFree(deviceInMat);

    return cudaStatus;
}

cudaError_t createIdentityMatrix(double* matrix, int dimension) {
    double* deviceMat = 0;
    cudaError_t cudaStatus = cudaSuccess;
    int matrixMemoryAllocationSize = pow(dimension, 2);

    cudaStatus = cudaMalloc((void**)&deviceMat, matrixMemoryAllocationSize * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        goto IdentityError;
    }

    cudaStatus = cudaMemcpy(deviceMat, matrix, matrixMemoryAllocationSize * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        goto IdentityError;
    }

    createIdentMat <<<dimension, dimension >> > (deviceMat, dimension);
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        goto IdentityError;
    }

    cudaStatus = cudaMemcpy(matrix, deviceMat, matrixMemoryAllocationSize * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        goto IdentityError;
    }


IdentityError:
    cudaFree(deviceMat);

    return cudaStatus;
}

cudaError_t multiplyMatrices(double* output, double* inputMatrix, int dimension) {
    double* deviceInputOne = 0;
    double* deviceInputTwo = 0;
    double* deviceOutput = 0;
    int matrixMemoryAllocationSize = pow(dimension, 2);
    cudaError_t cudaStatus = cudaSuccess;

    cudaStatus = cudaMalloc((void**)&deviceInputOne, matrixMemoryAllocationSize * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        goto MultiplyError;
    }

    cudaStatus = cudaMalloc((void**)&deviceInputTwo, matrixMemoryAllocationSize * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        goto MultiplyError;
    }

    cudaStatus = cudaMalloc((void**)&deviceOutput, matrixMemoryAllocationSize * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        goto MultiplyError;
    }

    cudaStatus = cudaMemcpy(deviceInputOne, output, matrixMemoryAllocationSize * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        goto MultiplyError;
    }

    cudaStatus = cudaMemcpy(deviceInputTwo, inputMatrix, matrixMemoryAllocationSize * sizeof(double), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        goto MultiplyError;
    }

    multiplyMatrix << <dimension, dimension >> > (deviceOutput, deviceInputOne, deviceInputTwo, dimension);
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        goto MultiplyError;
    }

    cudaStatus = cudaMemcpy(output, deviceOutput, matrixMemoryAllocationSize * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        goto MultiplyError;
    }

MultiplyError:
    cudaFree(deviceInputOne);
    cudaFree(deviceInputTwo);
    cudaFree(deviceOutput);

    return cudaStatus;
}

//Begin CUDA Function Implementations
__global__ void calcDistanceBetMats(double* outMat, double* inMatOne, double* inMatTwo, int dimension) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    //printf("calcDistanceBetMats [%d, %d]: %d\n", i, j, i + j * dimension);
    outMat[i + j * dimension] = std::pow(std::abs(inMatOne[i + j * dimension] - inMatTwo[i + j * dimension]),2);
    return;
}

__global__ void calcDistanceBetMatVec(double* outMat, double* inMat, double* inVec, int dimension) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = 0;
    //printf("calcDistanceBetMatVec [%d, %d]: %d\n", i, j, i + j * dimension);
    outMat[i + j * dimension] = std::pow(std::abs(inMat[i + j * dimension] - inVec[i + j * dimension]), 2);
    return;
}

__global__ void calcGaussCorr(double* outMat, double* inMat, int dimension, double variance, double a, double theta) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    //printf("calcGaussCorr [%d, %d]: %d\n", i, j, i + j * dimension);
    //Artifact of the cuComplex purge. Leaving here as it might be necessary later and will save a little time.
    double negOne = -1;
    outMat[i + j * dimension] = (variance - a) * std::exp(negOne * theta * inMat[i + j * dimension]);
    
    return;
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
            inputMatrix[x * n + y] /= inputMatrix[i * n + i];
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
                inputMatrix[x * n + y] -= inputMatrix[i * n + y] * inputMatrix[x * n + i];
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


//This might yield race condition errors
__global__ void multiplyMatrix(double* output, double* firInput, double* secInput, int dimension) {
    int x = blockIdx.x;
    int y = threadIdx.x;
    
    //Multiple each element pair and then sum them together
    for (int i = 0; i < dimension; i++) {
        output[x + y * dimension] += firInput[i + y * dimension] * secInput[x + i * dimension];
    }

    return;
}

__global__ void extendMat(double* outMat, double* inMat, int dimension) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int index = i + j * dimension;
    printf("extendMat [%d, %d]: %d\n", index, dimension, (dimension-1)*dimension );
    
    //If the current index is on the bottom-most row of the extended matrix
    if ( ((index + 1) % dimension) == 0) {
        printf("[%d] Hitting first statement\n", index);
        //If the sqrt(index+1) is equal to the dimension, it is the bottom right corner of the extended matrix
        if ((int)sqrt((float)(index + 1)) == dimension) {
            outMat[index] = 0.0;
        }
        //Otherwise, it is the bottom row of the extended matrix
        else {
            outMat[index] = 1.0;
        }
    }
    //If the current index is on the right-most column of the extended matrix. The previous statement should get the bottom right coordinate which is set to zero.
    else if ( index >= (dimension)*(dimension-1)) {
        printf("[%d] Hitting second statement\n", index);
        outMat[index] = 1.0;
    }
    //Otherwise, the index is within the original matrix
    else {
        printf("[%d] Hitting third statement\n", index);
        outMat[index] = inMat[index];
    }
    return;
}

__global__ void normalizeMatrix(double* outMat, double* inMat, double normalizingValue, int dimension) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    
    outMat[i + j * dimension] = inMat[i + j * dimension] / normalizingValue;

    return;
}

__global__ void createIdentMat(double* matrix, int dimension) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int index = i + j * dimension;

    if (index % (dimension + 1) == 0) {
        matrix[index] = 1;
    }
    else {
        matrix[index] = 0;
    }

    return;
}

void printMatrix(double inArray[], int dimension) {
    int index = 0;
    std::cout << "{";
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            index = i + (j * dimension);
            //printf("printMatrix[%d] = %f\n", index, inArray[index]);
            
            if (j + 1 < dimension) {
                printf("%0.60f,", inArray[index]);
            }
            else {
                printf("%0.60f", inArray[index]);
            }
            
        }
        
        if (i + 1 < dimension) {
            printf(";\n");
        }
        else {
            printf("}\n");
        }
        
    }
    return;
}