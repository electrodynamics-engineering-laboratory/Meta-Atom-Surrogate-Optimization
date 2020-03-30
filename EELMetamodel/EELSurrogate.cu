/*
File: EELSurrogate.cu
Date: 2-24-2020
Author: Joe Haun
Purpose: This file contains the EEL Surrogate class.  
*/

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>

//Function calculates the covariance of the given matrices
void calculateCovariance(double outputMatrix[], double firstInput[], double secondInput[], int dimension); 
double calculateAverage(double inputMatrix[], int dimension);

__global__ void covarianceHelper(double* outputMatrix, double* X, double meanX, double* Y, double meanY, int dimension){
	   int i = threadIdx.x;
	   int j = threadIdx.y;
       int index = i + (j*dimension);
	   outputMatrix[index] = (X[index]-meanX)*(Y[index]-meanY)/dimension;
	   return;
}

void calculateCovariance(double outputMatrix[], double firstInput[], double secondInput[], int dimension){
    //double firstInputMean = calculateAverage(firstInput, dimension);
    double array[4]  = { 1,2,3,4 };
    double firstInputMean = calculateAverage(array, dimension);
    double secondInputMean = calculateAverage(secondInput, dimension);
    covarianceHelper<<<dimension, dimension>>>(outputMatrix, firstInput, firstInputMean, secondInput, secondInputMean, dimension);
    return;
}

double calculateAverage(double inputMatrix[], int dimension){
    //Calculates the average value of the given matrix. Likely a way to do this with CUDA.
      double meanVal;
      double rowSum = 0;
      double colSum = 0;
      for(int i = 0; i < dimension; i++){
      	      colSum = 0;
      	      for(int j = 0; j < dimension; j++){
	      	      colSum += inputMatrix[i + (j*dimension)];	      
	      }
	      rowSum += colSum/dimension;
	      }
	meanVal = rowSum/dimension;
      return meanVal;
}