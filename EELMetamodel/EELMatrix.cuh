/*
File: EELMatrix.cuh
Date: 2-24-2020
Author: Joe Haun
Purpose: This file contains the EEL Matrix class definitions.  
*/

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>

class EELMatrix{
public:
	EELMatrix(int rows, int columns);
	~EELMatrix();
	EELMatrix* replicateMatrix(int rows, int columns);
	double at(int row, int column);
	void find(double value);

private:
	double* matrix;

}