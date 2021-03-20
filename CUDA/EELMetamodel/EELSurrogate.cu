/* @BEGIN_DOC_FILE!
File: EELSurrogate.cu
Purpose: This file contains the functions to implement a Kriging Surrogate model from an input data file.
Author: Joseph Haun
 @END_DOC_FILE!*/

//Surrogate Header File
#include "EELSurrogate.cuh"

//#include <cuComplex.h> //Unable to find proper documentation for cuComplex functionality. As such, some of the documentation relating to complex values is meaningless. 

//Begin Function Implementations
int colMajIndex(int x, int y, int dimension) {
    return x + y * dimension; #Dimension is columns
}

int rowMajIndex(int x, int y, int dimension) {
    return x * dimension + y; #Dimension is rows
}

double* readInputFile(std::string fileName) {
    //Assume no header lines
    return readInputFile(fileName, 0);
}

double* readInputFile(std::string fileName, int headerLines) {

    //Create filestream object and open target file
    std::fstream fileObject;
    fileObject.open(fileName, std::fstream::in);

    //Declare vector for file data
    std::vector<double> tempValues;

    //Declare vector and strings for parsing
    std::vector<std::string> row;
    std::string line, word, temp;

    //Declare counters
    int curRowIndex = 0;
    int curColIndex = 0;
    int numRows = 0;
    int numCols = 0;
    int headersSkipped = 0;

    while (fileObject >> line) {
        //Clear the vector of all values
        row.clear();

        //Break words apart
        if (headersSkipped >= headerLines) {
            std::stringstream inStream(line);
            
            while (getline(inStream, word, ',')) {
	    	//Attempt to convert the word to a double.
                try {
                    tempValues.push_back(std::stod(word.c_str()));
                }
                catch (const std::exception& e) {
                    tempValues.push_back(0.0);
                    printf("\nException occurred while attempting to convert '%s' to a double value. \nThe resulting output matrix has (%d, %d) set to zero.\n", word, curRowIndex, curColIndex);
                }
                
                curColIndex++;
                
            }

	    if(curColIndex > numCols){
	        numCols = curColIndex;
	    }	    
	    
            //Reset column index and increment row count
            numRows++;
            curColIndex = 0;
        }
        else {
            headersSkipped++;
        }
        
    }

    fileObject.close();

    //Allocate memory on heap for returned array and set all indices to 0.0
    double* fileValues = (double*)malloc(tempValues.size() * sizeof(double));
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            //printf("(%d,%d): %d -> %d\n", i, j, rowMajIndex(i,j,dimension), colMajIndex(i, j, fileLines));
            fileValues[colMajIndex(i, j, numRows)] = tempValues.at(rowMajIndex(i,j,numCols));
        }
    }

    return fileValues;
}

int* generateSample(int numRows, float size, bool random){

    //Generate the samples number for the array
    int numSamples = 0;
    //Enforce percentage limits
    if(size > 1.0 || size < 0.0){
       numSamples = numRows;
    }
    else{
       numSamples = numRows*size;
    }   

    //Allocate output memory
    int* samples (int*)malloc( numSamples* sizeof(int));
    
    //Initialize random seed
    std::srand(time(nullptr));

    bool duplicate = false;
    int checkCounter;
    for(int i = 0; i < numSamples; i++){
 	samples[i] = std::rand()%numRows;
	checkCounter = 0;
    	while(i > 0 && duplicate == false && checkCounter <= i){
	    //If the newly generated sample is already in the array
	    if(samples[checkCounter] == samples[i]){
	        //Generate new random value and reset counter
	        samples[i] = std::rand()%numRows;
		checkCounter = 0;
	    }
	    //Otherwise, increment counter
	    checkCounter++;
	}
    }
    return samples;
}

void parseRawData(double* outputParameters, double* outputData, double* inputData, int rows, int columns, int dataColumns, int* samples, int sampleCount){

    //If either pointer are not allocated, return
    if(outputParameters == nullptr || outputData == nullptr){
        return;
    }

    //For all sample rows
    for(int i = 0; i < sampleCount; i++){
        //For all columns in the row
        for(int j = 0; j < columns; j++){
            //If the column index is within the parameters column
            if(j < columns - dataColumns){
                //Copy the value to the parameters matrix
                outputParameters[colMajIndex(samples[i], j, columns-dataColumns)] = inputData[colMajIndex(samples[i], j, columns)];   
            }
            else{
                //Copy the value to the data matrix
                outputData[colMajIndex(samples[i], j, dataColumns)] = inputData[colMajIndex(samples[i], j, columns)];
            }
        }
    }
    return;
}

double kriging(std::string filename, int headerLines, int dataColumns, double theta, double variance, double nuggetEffect){
    //Begin variable definitions for data to be passed to GPU
    cudaError_t cudaStatus = cudaSuccess;
    int matrixMemoryAllocationSize = pow(dimension+1, 2);
    int vectorMemoryAllocationSize = dimension+1;
    double outputValue = 0;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess){
        goto KrigingError;
    }

    //Create a dynamic allocation of memory for the identity matrix and populate with values
    double* identityMatrix = (double*) malloc(matrixMemoryAllocationSize * sizeof(double));
    
    //Create a dynamic allocation of memory for a temporary holding matrix
    double* tempMatrixOne = (double*) malloc(extendedMatrixMemoryAllocationSize * sizeof(double));
    double* tempMatrixTwo = (double*) malloc(extendedMatrixMemoryAllocationSize * sizeof(double));
    for (int i = 0; i < pow(dimension+1,2); i++) {
        tempMatrixOne[i] = 0;
        tempMatrixTwo[i] = 0;
        identityMatrix[i] = 0;
    }

KrigingError:
   
    free(identityMatrix);
    free(tempMatrixOne);
    free(tempMatrixTwo);
    
    return outputValue;
}

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
        identityMatrix[i] = 0;
    }

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

    //Calculate distance between design sites and values at design sites. tempMatrixOne will hold the output
    cudaStatus = calculateDistanceBetweenMatrices(tempMatrixOne, designSite, designSiteValues, dimension);
    if (cudaStatus != cudaSuccess) {
        goto SetupError;
    }

    //Calculate distance between test site and design site. tempMatrixTwo will hold the output
    cudaStatus = calculateDistanceBetweenMatrixVector(tempMatrixTwo, designSite, testSite, dimension);
    if (cudaStatus != cudaSuccess) {
        goto SetupError;
    }

    //Calculate the covariance between design sites and values at design sites. tempMatrixOne will hold the output
    cudaStatus = calculateGaussianCorrelation(tempMatrixOne, tempMatrixOne, variance, a, theta, dimension);
    if (cudaStatus != cudaSuccess) {
        goto SetupError;
    }

    //Calculate the covariance between test sites and design sites. tempMatrixTwo will hold the output
    cudaStatus = calculateGaussianCorrelation(tempMatrixTwo, tempMatrixTwo, variance, a, theta, dimension);
    if (cudaStatus != cudaSuccess) {
        goto SetupError;
    }

    //Clear the extraneous elements of the matrix after the valid vector elements to prevent later calculation errors.
    for (int i = dimension; i < pow(dimension + 1, 2); i++) {
        tempMatrixTwo[i] = 0;
    }

    //Extend the covariance matrix between the design site and design site values
    cudaStatus = extendMatrix(tempMatrixOne, tempMatrixOne, dimension);
    if (cudaStatus != cudaSuccess) {
        goto SetupError;
    }

    //Extend the covariance vector between test site and design sites. 
    tempMatrixTwo[dimension + 0*dimension] = 1; //Add 1 to the last row of the matrix, unclear if this is correct
    
    //Calculate inverse of extended covariance matrix
    cudaStatus = invertMatrix(identityMatrix, tempMatrixOne, dimension + 1);
    if (cudaStatus != cudaSuccess) {
        goto SetupError;
    }

    //Calculate extended weights vector, tempMatrixTwo will hold the result
    cudaStatus = calculateWeightVector(tempMatrixTwo, tempMatrixOne, tempMatrixTwo, dimension + 1);
    if (cudaStatus != cudaSuccess) {
        goto SetupError;
    }

    //Calculate estimate value at test site, the ultimate output. 
    cudaStatus = multiplyMatrices(tempMatrixTwo, designSiteValues, dimension + 1); //Only consider elements within the dimension as the final value in the weight matrix is the lamda value (AKA not needed)
    if (cudaStatus != cudaSuccess) {
        goto SetupError;
    }

    //Grab the first, and only, value of tempMatrixTwo as the ultimate output value
    outputValue = tempMatrixTwo[0];

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
    cudaStatus = cudaMemcpy(outputMatrix, deviceOutMat, matrixMemoryAllocationSize * sizeof(double), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        goto ExtendError;
    }

ExtendError:
    cudaFree(deviceOutMat);
    cudaFree(deviceInMat);

    return cudaStatus;
}

cudaError_t invertMatrix(double* outputMatrix, double* inputMatrix, int dimension) {
    
    double* deviceInMat = 0;
    double* deviceIdenMat = 0;
    double* deviceBuffer = 0;
    bool* deviceFlag = 0;
    bool* flag = new bool;
    *flag = true;

    int matrixMemoryAllocationSize = pow(dimension, 2);
    cudaError_t cudaStatus = cudaSuccess;

    cudaStatus = cudaMalloc((void**)&deviceBuffer, 2*dimension * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        goto InvertError;
    }

    cudaStatus = cudaMalloc((void**)&deviceFlag, sizeof(bool));
    if (cudaStatus != cudaSuccess) {
        goto InvertError;
    }

    cudaStatus = cudaMalloc((void**)&deviceInMat, matrixMemoryAllocationSize * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        goto InvertError;
    }

    cudaStatus = cudaMalloc((void**)&deviceIdenMat, matrixMemoryAllocationSize * sizeof(double));
    if (cudaStatus != cudaSuccess) {
        goto InvertError;
    }

    cudaStatus = cudaMemcpy(deviceFlag, flag, sizeof(bool), cudaMemcpyHostToDevice);
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
        resetBuffers << <dimension, dimension >> > (deviceBuffer, deviceFlag, dimension);
        normalizeRows << <dimension, dimension >> > (deviceIdenMat, deviceInMat, deviceBuffer, deviceFlag, dimension, i);
        pivotDown << <dimension, dimension >> > (deviceIdenMat, deviceInMat, dimension, i);
        
    }

    for (int i = dimension - 1; i > 0; i--) {
        resetBuffers << <dimension, dimension >> > (deviceBuffer, deviceFlag, dimension);
        pivotUp << <dimension, dimension >> > (deviceIdenMat, deviceInMat, deviceBuffer, deviceFlag, dimension, i);
    }

    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        goto InvertError;
    }

    //Synchronize device to enable a copy of the result
    cudaStatus = cudaDeviceSynchronize();
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
    cudaFree(deviceBuffer);
    cudaFree(deviceFlag);
    free(flag);

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
    cudaStatus = cudaMemcpy(deviceInvCovMat, invertedCovarianceMatrix, matrixMemoryAllocationSize * sizeof(double), cudaMemcpyHostToDevice);
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

void printMatrix(double inArray[], int dimension) {
    int index = 0;
    std::cout << "{";
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            index = colMajIndex(i,j,dimension);
            //printf("printMatrix[%d] = %f\n", index, inArray[index]);

            if (j + 1 < dimension) {
                printf("%0.5f,", inArray[index]);
            }
            else {
                printf("%0.5f", inArray[index]);
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

void printMatrix(double inArray[], int numRows, int numColumns) {
    
    int index = 0;
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numColumns; j++) {
            index = colMajIndex(i,j,numRows);
            //printf("printMatrix[%d] = %f\n", index, inArray[index]);

            if (j + 1 < numColumns) {
                printf("%0.5f,", inArray[index]);
            }
            else {
                printf("%0.5f", inArray[index]);
            }

        }

        if (i + 1 < numRows) {
            printf(";\n");
        }
        else {
            printf("}\n");
        }

    }
    return;

}

//Begin CUDA Function Implementations
__global__ void calcDistance(double* outMat, double* inMat, int inRows, int inCols) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    //Generate thread index based on input columns
    int index = i + j * inCols;
    //Calculate specific i and j parameters based on input rows/cols
    i = index % inCols;
    j = index / inCols;
    //For all columns in the input matrix, 
    for(int k = 0; k < inRows){
    	//Add the difference of the values to the output matrix
        outMat[index] += std::pow(inMat[k + i*inCols] - inMat[k + j*inCols], 2);
    }
    return;
}

__global__ void calcGaussCorr(double* outMat, double* inMat, int dimension, double variance, double a, double theta) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int index = i + j * dimension;
    //Artifact of the cuComplex purge. Leaving here as it might be necessary later and will save a little time.
    double negOne = -1;
    outMat[index] = (variance - a) * std::exp(negOne * theta * inMat[index]);
    
    return;
}

__global__ void resetBuffers(double* vals, bool* flag, int dimension) {
    
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int index = i + j * dimension;

    //Reset the flag to false if the zeroth index. This prevents race conditions between threads as there is only one value to write. 
    if (index == 0) {
        *flag = false;
    }

    //If the index is below double the dimensions, set the value equal to zero. Not strictly necessary, but prevents data manipulation by extra threads in the case of thread optimization.
    if (index < 2*dimension) {
        vals[index] = 0;
    }
    return;
}

__global__ void normalizeRows(double* idenMat, double* inMat, double* firstVals, bool* flag, int dimension, int targetCol) {
    /*BEGIN DOC
    Function:       normalizeRows(idenMat, inMat, dimension, flag, targetCol, firstVals)
    Description:    Perform the necessary normalization for the pivot down step of the Gauss-Jordan Elimination algorithm.
    Inputs:         identMat (double*) - A pointer to an array of double values that represent an identity matrix. Stored in row major format
                    inMat (double*) - A pointer to an array of double values that represent a matrix of values to invert. Stored in row major format.
                    dimension (int) - The dimension, n, of the square matrices.
                    flag (bool*) - A pointer to a boolean that functions as the flag to signal all threads to complete operations.
                    targetCol (int) - The target column index from which the pivot downwards operation occurs.
    Outputs:        None
    Notes:          Each row will be normalized by the value in the target column. The target column value is also used as a row index modifier to prevent the further normalization of values on the diagonal that have a lower index that the leading value of the targeted column.
    END DOC*/

    //Get the index value for the particular thread by 
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int index = i + j * dimension;
    
    //If the thread index is that of the first row and targeted column
    if (index == targetCol * (dimension + 1)) {
        //Get all subsequent element in the column of index
        for (int k = index % dimension; k < dimension; k++) {
            firstVals[k] = inMat[index - targetCol + k];
        }
        *flag = true;
    }

    //Use the flag as a lock to prevent divide-by-zero errors
    while (*flag == false) {
        //Simply wait for the array of first values to be set. 
    }

    //Grab the normalizing value for the particular thread by dividing the index by dimension in integer division. Decimal values should not be present.
    double normVal = firstVals[index % dimension];
    //If normalizing value is not zero and the index falls within the desired submatrix
    if (normVal != 0 && index >= targetCol % dimension) { //normVal is a double, might not ever equal 0. Might need to figure out a better way to exclude values close to zero. Establish precision limits.
        //Divide the identity and input matrix position by the normalizing value
        idenMat[index] /= normVal;
        inMat[index] /= normVal;
    }

    return;
}
    
__global__ void pivotDown(double* idenMat, double* inMat, int dimension, int targetRow) {
    /*BEGIN DOC
    Function:       pivotDown(idenMat, inMat, dimension, flag, targetRow)
    Description:    Perform the pivot down step of the Gauss-Jordan Elimination algorithm. 
    Inputs:         identMat (double*) - A pointer to an array of double values that represent an identity matrix. Stored in row major format
                    inMat (double*) - A pointer to an array of double values that represent a matrix of values to invert. Stored in row major format.
                    dimension (int) - The dimension, n, of the square matrices. 
                    targetRow (int) - The target row index from which the pivot downwards operation occurs.
    Outputs:        None
    Notes:          Every row below the target row will have their column values subtracted by the pivot row's respective column value. The pivot row remains the same. 
    END DOC*/

    //Get the index value for the particular thread by multiplying the thread index values with the dimension in a row major format
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int index = i + j * dimension;

    //If the current index falls above the targeted row and the leading element is not zero, subtract the target row column value from the index's value
    if ( (index % dimension > targetRow) && (inMat[index % dimension + targetRow*dimension] != 0)) {
        //double temp = inMat[index - (index % dimension - targetRow)];
        int temp = index - (index % dimension - targetRow);
        idenMat[index] -= idenMat[temp];
        inMat[index] -= inMat[temp];
    }
    return;
}

__global__ void pivotUp(double* idenMat, double* inMat, double* lastVals, bool* flag, int dimension, int targetRow) {
    /*BEGIN DOC
    Function:       pivotUp(idenMat, inMat, lastVals, flag, dimension, targetRow)
    Description:    Perform the pivot up step of the Gauss-Jordan Elimination algorithm.
    Inputs:         identMat (double*) - A pointer to an array of double values that represent an identity matrix. Stored in row major format
                    inMat (double*) - A pointer to an array of double values that represent a matrix of values to invert. Stored in row major format.
                    dimension (int) - The dimension, n, of the square matrices.
                    flag (bool*) - A pointer to a boolean that functions as the flag to signal all threads to complete operations.
                    targetRow (int) - The target row index from which the pivot upwards operation occurs. 
                    lastVals (double*) - A pointer to an array that will be used to store values for operations.
    Outputs:        None
    Notes:          Every row above the target row will have their column values subtracted by the pivot row's respective column value multiplied by the column value indicated by the targetRow value for that row. As the matrices are square, the target row is also the target column. The pivot row remains the same.
    END DOC*/

    //Get the index value for the particular thread by 
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int index = i + j * dimension;

    //From the leading pivot column, find non-zero values from this row to the zeroth row
    if (index == targetRow * (dimension + 1)) {
        for (int k = index % dimension; k >= 0; k--) {
            lastVals[(index - k)%dimension] = inMat[index - k];
        }
        *flag = true;
    }

    while (*flag == false) {
        //Wait for the lastVals array to be populated
    }

    //If the index falls in a row below the target, perform subtraction operations
    if (index % dimension < targetRow ) {
        //Multiply the index value by the trailing row value and subtract value from current index in both matrices
        idenMat[index] -= lastVals[index % dimension] * idenMat[index + targetRow - (index % dimension)];
        inMat[index] -= lastVals[index % dimension] * inMat[index + targetRow - (index % dimension)];
    }
    return;
}

//This might yield race condition errors
__global__ void multiplyMatrix(double* output, double* firInput, double* secInput, int dimension) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    //Generate index for the output based on thread and block parameters
    int index = i + j * dimension;
    
    //Calculate target row and column values from the index to allow for dot product to occur
    int row = index % dimension;
    int col = index / dimension;

    //Multiply each element pair and then sum them together
    for (int k = 0; k < dimension; k++) {
        output[index] += firInput[row + k * dimension] * secInput[k +  col * dimension];
    }

    return;
}

__global__ void extendMat(double* outMat, double* inMat, int dimension) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int index = i + j * dimension;
    
    //If the current index is on the bottom-most row of the extended matrix
    if ( ((index + 1) % dimension) == 0) {
        
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
        
        outMat[index] = 1.0;
    }
    //Otherwise, the index is within the original matrix
    else {
        
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