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
void calculateCovariance(double* outputMatrix, double* firstInput, double* secondInput, int dimension); 
double calculateAverage(double inputMatrix, int dimension);

__global__ void covarianceHelper(double* outputMatrix, double* X, double meanX, double* Y, double meanY, int dimension){
	   int i = threadIdx.x;
	   int j = threadIdx.y;
	   outputMatrix[i][j] = (X[i][j]-meanX)*(Y[i][j]-meanY)/dimension;
	   return;
}

void calculateCovariance(double* outputMatrix, double* firstInput, double* secondInput, int dimension){
     double firstInputMean = calculateAverage(firstInput, dimension);
     double secondInputMean = calculateAverage(secondInput, dimension);
     covarianceHelper<<<dimension, dimension>>>(outputMatrix, firstInput, firstInputMean, secondInput, secondInputMean, dimension);
     
}

double calculateAverage(double inputMatrix, int dimension){
      double meanVal;
      double rowSum = 0;
      double colSum = 0;
      for(int i = 0; i < dimension; i++){
      	      colSum = 0;
      	      for(int j = 0; j < dimension; j++){
	      	      colSum += inputMatrix[i][j];	      
	      }
	      rowSum += colSum/sizeY;
	      }
	meanVal = rowSum/sizeX;
      return meanVal;
}

int main(){

    return 0;
}