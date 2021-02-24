/* File:	EELSurrogate.cuh
   Purpose:	Header file for EELSurrogate.cu
   Created:	2-7-2021
   Modified:	2-7-2021
   Notes: 
*/

//Standard C++ libraries
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <math.h>
//#include <conio.h>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

#ifndef SURROGATE_H
#define SURROGATE_H

//Needed to provide full path for CUDA functions on a Windows System. This should be changed for the local machine, unless Visual Studio is setup correctly.
//CUDA libraries
//#include "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.0\include\cuda_runtime.h"
//#include "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.0\include\device_launch_parameters.h"
//#include <C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.0\include\cuda_fp16.h>

//Begin CPU Function Definitions
/* @BEGIN_DOC_FUNC!
Function:    colMajIndex(int x, int y, int dimension)
Purpose:     Provide a CPU version of a column-major index given two-dimensional coordinates.
Inputs:      x (int) - The x coordinate of the matrix
             y (int) - The y coordinate of the matrix
             dimension (int) - The dimensionality of the matrix. 
Outputs:     coordinate (int) - The coordinate in a flat column-major matrix.
@END_DOC_FUNC! */
int colMajIndex(int x, int y, int dimension);

/* @BEGIN_DOC_FUNC!
Function:    rowMajIndex(int x, int y, int dimension)
Purpose:     Provide a CPU version of a row-major index given two-dimensional coordinates.
Inputs:      x (int) - The x coordinate of the matrix
             y (int) - The y coordinate of the matrix
             dimension (int) - The dimensionality of the matrix.
Outputs:     coordinate (int) - The coordinate in a flat column-major matrix.
@END_DOC_FUNC! */
int rowMajIndex(int x, int y, int dimension);

/* @BEGIN_DOC_FUNC!
Function:    readInputFile(std::string fileName)
Purpose:     Read a CSV file and format the data into a column-major formatted matrix. Automatically assumes zero header lines to skip.
Inputs:      fileName (std::string) - The file path of the CSV to open and read.
Outputs:     fileValues (double*) - The values within a flat column-major formatted matrix.
Notes:       The returned matrix is in column-major format.
@END_DOC_FUNC! */
double* readInputFile(std::string fileName);

/* @BEGIN_DOC_FUNC!
Function:    readInputFile(std::string fileName, int headerLines)
Purpose:     Read a CSV file and format the data into a column-major formatted matrix. 
Inputs:      fileName (std::string) - The file path of the CSV to open and read.
             headerLines (int) - An integer that indicates the number of header lines to skip. 
Outputs:     fileValues (double*) - The estimator value that is calculated for the given target, design sites, and design site values.
Notes:       The returned matrix is in column-major format.
@END_DOC_FUNC! */
double* readInputFile(std::string fileName, int headerLines);

/* @BEGIN_DOC_FUNC!
Function:    generateSample(int numRows, float size, bool random)
Purpose:     Create an array of sample indices to be chosen from the main data matrix.
Inputs:      numRows (int) - The number of rows from which to choose.
	     size (float) - The percentage of the total rows to use for a sample. 
	     random (bool) - Determine if the sample is random or not.
Outputs:     sampleChoices (std::vector<int>) - A vector of indices to use for the sample.
Notes:       The returned matrix is in column-major format.
@END_DOC_FUNC! */
std::vector<int> generateSample(int numRows, float size, bool random);

/* @BEGIN_DOC_FUNC!
Function:    appendRow(double* inputMatrix, double* outputMatrix, int inputRow, int outputRow)
Purpose:     Copy one row to another row.
Inputs:	     inputMatrix (double*) - The matrix to copy the row from.
	     outputMatrix (double*) - The matrix to append the row.
	     inputRow (int) - The target row from the input matrix.
	     outputRow (int) - The target row in the output matrix.
Outputs:     sampleChoices (int*) - An array of indices to use for the sample.
Notes:       The returned matrix is in column-major format.
@END_DOC_FUNC! */
void appendRow(double* inputMatrix, double* outputMatrix, int inputRow, int outputRow);

/* @BEGIN_DOC_FUNC!
Function:    metamodelSetup(int dimension, double theta, double variance, double a, double* designSite, double* testSite, double* designSiteValues)
Purpose:     Set up and perform calculations on the GPU that relate to a Kriging metamodel
Inputs:      dimension (int) - An integer value representing the number of rows for all matrices and vectors and the number columns for all matrices
             theta (double) - (THIS MIGHT BE AN INCORRECT DESCRIPTION) A value in radians that represents the angle between the two matrices, this value should only ever be real
             variance (double) - The value of the variance to use in calculations.
             a (double) - The nugget effect value to use in calculations.
             designSite (double*) - a dimension-by-dimension matrix that represents the datapoints previously collected to be used to estimate a new value based on test site
             testSite (double*)- a dimension-by-dimension matrix with only the first column that make up the vector that represents the points at which a new value is to be calculated
             designSiteValues (double*) - a dimension-by-dimension matrix that represent the values at particular design sites. Each matrix has a 1-to-1 correlation of site to value (i.e. designSite[0] is paired with designSiteValues[0])
Outputs:     outputValue (double) - The estimator value that is calculated for the given target, design sites, and design site values.
Notes:       All matrices are in column-major format.
@END_DOC_FUNC! */
//Function sets up and performs calculations on the GPU
double metamodelSetup(int dimension, double theta, double variance, double a, double* designSite, double* testSite, double* designSiteValues );

/* @BEGIN_DOC_FUNC!
Function:    calculateGaussianCorrelation(double* outputMatrix, double* inMatrix, double variance, double a, double theta, int dimension)
Purpose:     Set up the GPU and then calculate covariance of two matrices through Gaussian Correlation.
Inputs:      outputMatrix (double*) - a dimension-by-dimension matrix. This input will be overwritten with the calculated matrix. 
             inMatrix (double*)- a dimension-by-dimension matrix of values with which to calculate the Gaussian Correlation.
             theta (double) - (THIS MIGHT BE AN INCORRECT DESCRIPTION) a value in radians that represents the angle between the two matrices, this value should only ever be real.
             variance (double) - the value of the variance for use in calculation of the Gaussian Correlation.
             dimension (int) - an integer value representing the number of rows for all matrices and vectors and the number columns for all matrices
Outputs:     cudaStatus (cudaError_t) - The error status of the CUDA operations.
Notes:       All matrices are in column-major format.
@END_DOC_FUNC! */
cudaError_t calculateGaussianCorrelation(double* outputMatrix, double* inMatrix, double variance, double a, double theta, int dimension);

/* @BEGIN_DOC_FUNC!
Function:    calculateDistanceBetweenMatrices(outputMatrix, inMatrixOne, inMatrixTwo, dimension)
Purpose:     Set up the GPU and then calculate the distance between two matrices.
Inputs:      outputMatrix (double*)- a dimension-by-dimension matrix. This input will be overwritten with the calculated matrix.
             inMatrixOne (double*) - a dimension-by-dimension matrix in column-major format.
             inMatrixTwo (double*) - a dimension-by-dimension matrix in column-major format. 
             dimension (int) - an integer value representing the number of rows for all matrices and vectors and the number columns for all matrices
Outputs:     cudaStatus (cudaError_t) - the error status of the CUDA operations.
Notes:       All matrices are in column-major format.
@END_DOC_FUNC! */
cudaError_t calculateDistanceBetweenMatrices(double* outputMatrix, double* inMatrixOne, double* inMatrixTwo, int dimension);

/* @BEGIN_DOC_FUNC!
Function:    calculateDistanceBetweenMatrixVector(double* outputMatrix, double* inMatrix, double* inVector, int dimension)
Purpose:     Set up the GPU and then calculate the distance between a matrix and a vector.
Inputs:      outputMatrix (double*)- a dimension-by-dimension matrix. This input will be overwritten with the calculated matrix.
             inMatrix (double*) - a dimension-by-dimension matrix in column-major format.
             inVector (double*) - a dimension-by-dimension matrix in column-major format with only the first column filled with values to represent the vector. 
             dimension (int) - an integer value representing the number of valid rows for all matrices and vectors and the number of valid columns for all matrices
Outputs:     cudaStatus (cudaError_t) - The error status of the CUDA operations.
Notes:       All matrices are in column-major format. While the inputs for vectors have the same dimensionality as matrices, always ensure that the input vectors have all values outside 0 to dimension-1 set to zero. These indices are not ignored in calculations, but will not affect the calculations when set to zero.
@END_DOC_FUNC! */
cudaError_t calculateDistanceBetweenMatrixVector(double* outputMatrix, double* inMatrix, double* inVector, int dimension);

//Function sets up GPU for calculating the weights vector
/* @BEGIN_DOC_FUNC!
Function:    calculateWeightVector(double* outputVector, double* invertedCovarianceMatrix, double* covarianceVector, int dimension)
Purpose:     Function sets up and performs calculations on the GPU that relate to a metamodel
Inputs:      outputVector (double*)- a dimension-by-dimension matrix. This input will be overwritten with the calculated vector.
             invertedCovarianceMatrix (double*) - a dimension-by-dimension matrix in column-major format.
             covarianceVector (double*) - a dimension-by-dimension matrix in column-major format with only the first column filled with values to represent the vector. 
             dimension (int) - an integer value representing the number of valid rows for all matrices and vectors and the number of valid columns for all matrices
Outputs:     cudaStatus (cudaError_t) - The error status of the CUDA operations.
Notes:       All matrices are in column-major format. While the inputs for vectors have the same dimensionality as matrices, always ensure that the input vectors have all values outside 0 to dimension-1 set to zero. These indices are not ignored in calculations, but will not affect the calculations when set to zero.
@END_DOC_FUNC! */
cudaError_t calculateWeightVector(double* outputVector, double* invertedCovarianceMatrix, double* covarianceVector, int dimension);


/* @BEGIN_DOC_FUNC!
Function:    extendMatrix(double* outputMatrix, double* inputMatrix, int dimension)
Purpose:     Sets up GPU and then extends a matrix
Inputs:      outputMatrix (double*)- a (dimension+1)-by-(dimension+1) matrix. This input will be overwritten with the calculated matrix.
             inputMatrix (double*) - a dimension-by-dimension matrix in column-major format.
             dimension (int) - an integer value representing the number of valid rows for all matrices and vectors and the number of valid columns for all matrices
Outputs:     cudaStatus (cudaError_t) - The error status of the CUDA operations.
Notes:       All matrices are in column-major format.
@END_DOC_FUNC! */
cudaError_t extendMatrix(double* outputMatrix, double* inputMatrix, int dimension);

//Function sets up GPU for calculating the inverse of a matrix using Gauss-Jordan elimination
/* @BEGIN_DOC_FUNC!
Function:    invertMatrix(double* outputMatrix, double* inputMatrix, int dimension)
Purpose:     Function sets up and performs calculations on the GPU that relate to a metamodel
Inputs:      outputMatrix (double*)- a dimension-by-dimension identity matrix in column-major format. This input will be overwritten with the calculated matrix.
             inputMatrix (double*) - a dimension-by-dimension matrix in column-major format.
             dimension (int) - an integer value representing the number of valid rows for all matrices and vectors and the number of valid columns for all matrices
Outputs:     cudaStatus (cudaError_t) - The error status of the CUDA operations.
Notes:       All matrices are in column-major format. It is required of this function that the output matrix be an identity matrix. Otherwise, the Gauss-Jordan elimination algorithm is meaningless. 
@END_DOC_FUNC! */
cudaError_t invertMatrix(double* outputMatrix, double* inputMatrix, int dimension);

/* @BEGIN_DOC_FUNC!
Function:    createIdentityMatrix(double* matrix, int dimension)
Purpose:     Sets the given matrix to the identity matrix.
Inputs:      matrix (double*) - A dimension-by-dimension matrix.
             dimension (int) - an integer value representing the number of valid rows and columns for the matrix.
Outputs:     cudaStatus (cudaError_t) - The error status of the CUDA operations.
Notes:       All matrices are in column-major format.
@END_DOC_FUNC! */
cudaError_t createIdentityMatrix(double* matrix, int dimension);

/* @BEGIN_DOC_FUNC!
Function:    multiplyMatrices(double* output, double* inputMatrix, int dimension)
Purpose:     Set up the GPU and then multiply two matrices.
Inputs:      output (double*) - a dimension-by-dimension matrix to be used as one input in the multiplication. This matrix will be overwritten with the calculated matrix.
             inputMatrix (double*) - a dimension-by-dimension matrix to be used as on input in the multiplication.
             dimension (int) - an integer value representing the number of valid rows for all matrices and vectors and the number of valid columns for all matrices
Outputs:     cudaStatus (cudaError_t) - The error status of the CUDA operations.
Notes:       All matrices are in column-major format.
@END_DOC_FUNC! */
cudaError_t multiplyMatrices(double* output, double* inputMatrix, int dimension);

//Testing function to print out arrays after computation.
/* @BEGIN_DOC_FUNC!
Function:    printMatrix(double inArray[], int dimension)
Purpose:     Function sets up and performs calculations on the GPU that relate to a metamodel
Inputs:      inArray (double[]) - The array of values that will be printed out as a matrix.
             dimension (int) - An integer value representing the number of rows for all matrices and vectors and the number columns for all matrices.
Outputs:     None
Notes:       All matrices are in column-major format.
@END_DOC_FUNC! */
void printMatrix(double inArray[], int dimension);

void printMatrix(double inArray[], int numRows, int numColumns);

//Begin GPU Function Definitions
/* @BEGIN_DOC_FUNC!
Function:    calcDistanceBetMats(double* outMat, double* inMatOne, double* inMatTwo, int dimension)
Purpose:     Calculate distance between two matrices on the GPU.
Inputs:      outMat (double*) - a dimension-by-dimension matrix to be used as the output for calculations.
             inMatOne (double*) - a dimension-by-dimension matrix.
             inMatTwo (double*) - a dimension-by-dimension matrix.
             dimension (int) - an integer value representing the number of valid rows for all matrices and vectors and the number of valid columns for all matrices
Outputs:     None
Notes:       All matrices are in column-major format.
@END_DOC_FUNC! */
__global__ void calcDistanceBetMats(double* outMat, double* inMatOne, double* inMatTwo, int dimension);

/* @BEGIN_DOC_FUNC!
Function:    calcDistanceBetMatVec(double* outMat, double* inMat, double* inVec, int dimension)
Purpose:     Calculate the distance between a matrix and a vector.
Inputs:      outMat (double*) - a dimension-by-dimension matrix to be used as the output for calculations.
             inMat (double*) - a dimension-by-dimension matrix.
             inVec (double*) - a dimension-by-dimension matrix with only the first column filled with active values.
             dimension (int) - an integer value representing the number of valid rows for all matrices and vectors and the number of valid columns for all matrices
Outputs:     None
Notes:       All matrices are in column-major format.
@END_DOC_FUNC! */
__global__ void calcDistanceBetMatVec(double* outMat, double* inMat, double* inVec, int dimension);

/* @BEGIN_DOC_FUNC!
Function:    calcGaussCorr(double* outMat, double* inMat, int dimension, double variance, double a, double theta)
Purpose:     Calculate the Gaussian Correlation of two matrices.
Inputs:      outMat (double*) - A dimension-by-dimension matrix to be used as the output for calculations.
             inMat (double*) - A dimension-by-dimension matrix.
             variance (double) - The variance value to be used in calculating the Gaussian Correlation.
             a (double) - The nugget effect value to be used in calculations.
             theta (double) - The theta value to be used in calculations.
             dimension (int) - an integer value representing the number of rows for all matrices and vectors and the number columns for all matrices
Outputs:     None
Notes:       All matrices are in column-major format.
@END_DOC_FUNC! */
__global__ void calcGaussCorr(double* outMat, double* inMat, int dimension, double variance, double a, double theta);

/* @BEGIN_DOC_FUNC!
Function:    normalizeMatrix(double* outMat, double* inMat, double normalizingValue, int dimension)
Purpose:     Normalize a given matrix to a value.
Inputs:      outMat (double*) - A dimension-by-dimension matrix to be used as the output for calculations.
             inMat (double*) - A dimension-by-dimension matrix.
             normalizingValue (double) - The value with which to normalize the input matrix.
             dimension (int) - an integer value representing the number of rows for all matrices and vectors and the number columns for all matrices
Outputs:     None
Notes:       All matrices are in column-major format.
@END_DOC_FUNC! */
__global__ void normalizeMatrix(double* outMat, double* inMat, double normalizingValue, int dimension);

//CUDA function to extend a given matrix
/* @BEGIN_DOC_FUNC!
Function:    extendMat(double* outMat, double* inMat, int dimension)
Purpose:     Function sets up and performs calculations on the GPU that relate to a metamodel
Inputs:      outMat (double*) - A dimension-by-dimension matrix to be used as the output for calculations.
             inMat (double*) - A dimension-by-dimension matrix.
             dimension (int) - an integer value representing the number of rows for all matrices and vectors and the number columns for all matrices
Outputs:     None
Notes:       All matrices are in column-major format.
@END_DOC_FUNC! */
__global__ void extendMat(double* outMat, double* inMat, int dimension);

//CUDA functions to calculate the inverse of a given matrix using Gauss-Jordan elimination method
/* @BEGIN_DOC_FUNC!
Function:    resetBuffers(double* vals, bool* flag, int dimension)
Purpose:     Function sets up and performs calculations on the GPU that relate to a metamodel
Inputs:      vals (double*) - A dimension-by-1 matrix buffer to be used to store values for calculations.
             flag (bool*) - A boolean flag to be used to determine calculation order.
             dimension (int) - an integer value representing the number of rows for all matrices and vectors and the number columns for all matrices
Outputs:     None
Notes:       All matrices are in column-major format.
@END_DOC_FUNC! */
__global__ void resetBuffers(double* vals, bool* flag, int dimension);

/* @BEGIN_DOC_FUNC!
Function:    normalizeRows(double* idenMat, double* inMat, double* firstVals, bool* flag, int dimension, int targetCol)
Purpose:     Function sets up and performs calculations on the GPU that relate to a metamodel
Inputs:      idenMat (double*) - A dimension-by-dimension matrix to be used as the output for calculations.
             inMat (double*) - A dimension-by-dimension matrix.
             firstVals (double*) - A buffer array to be used to prevent calculation errors.
             flag (bool*) - A flag to prevent calculations until firstVals is populated.
             dimension (int) - An integer value representing the number of rows for all matrices and vectors and the number columns for all matrices.
             targetCol (int) - The column to target for the normalize operations.
Outputs:     None
Notes:       All matrices are in column-major format.
@END_DOC_FUNC! */
__global__ void normalizeRows(double* idenMat, double* inMat, double* firstVals, bool* flag, int dimension, int targetCol);

/* @BEGIN_DOC_FUNC!
Function:    pivotDown(double* idenMat, double* inMat, int dimension, int targetRow)
Purpose:     Function sets up and performs calculations on the GPU that relate to a metamodel
Inputs:      idenMat (double*) - A dimension-by-dimension matrix to be used as the output for calculations.
             inMat (double*) - A dimension-by-dimension matrix.
             dimension (int) - An integer value representing the number of rows for all matrices and vectors and the number columns for all matrices.
             targetRow (int) - The row to target for the pivot down operations.
Outputs:     None
Notes:       All matrices are in column-major format.
@END_DOC_FUNC! */
__global__ void pivotDown(double* idenMat, double* inMat, int dimension, int targetRow);

/* @BEGIN_DOC_FUNC!
Function:    pivotUp(double* idenMat, double* inMat, double* lastVals, bool* flag, int dimension, int targetRow)
Purpose:     Function sets up and performs calculations on the GPU that relate to a metamodel
Inputs:      idenMat (double*) - A dimension-by-dimension matrix to be used as the output for calculations.
             inMat (double*) - A dimension-by-dimension matrix.
             lastVals (double*) - A buffer array to be used to prevent calculation errors.
             flag (bool*) - A flag to prevent calculations until firstVals is populated.
             dimension (int) - An integer value representing the number of rows for all matrices and vectors and the number columns for all matrices.
             targetRow (int) - The row to target for the pivot down operations.
Outputs:     None
Notes:       All matrices are in column-major format.
@END_DOC_FUNC! */
__global__ void pivotUp(double* idenMat, double* inMat, double* lastVals, bool* flag, int dimension, int targetRow);

//CUDA function to calculate the dot product of two inputs
/* @BEGIN_DOC_FUNC!
Function:    multiplyMatrix(double* output, double* firInput, double* secInput, int dimension)
Purpose:     Function sets up and performs calculations on the GPU that relate to a metamodel
Inputs:      output (double*) - A dimension-by-dimension matrix to be used as the output for calculations.
             firInput (double*) - A dimension-by-dimension matrix.
             secInput (double*) - A dimension-by-dimension matrix.
             dimension (int) - an integer value representing the number of rows for all matrices and vectors and the number columns for all matrices
Outputs:     None
Notes:       All matrices are in column-major format.
@END_DOC_FUNC! */
__global__ void multiplyMatrix(double* output, double* firInput, double* secInput, int dimension);

//CUDA function to create an identity matrix of the given dimension
/* @BEGIN_DOC_FUNC!
Function:    createIdentMat(double* matrix, int dimension)
Purpose:     Function sets up and performs calculations on the GPU that relate to a metamodel
Inputs:      matrix (double*) - A matrix to be overwritten with the identity matrix.
             dimension (int) - an integer value representing the number of rows for all matrices and vectors and the number columns for all matrices
Outputs:     None
Notes:       All matrices are in column-major format.
@END_DOC_FUNC! */
__global__ void createIdentMat(double* matrix, int dimension);

#endif