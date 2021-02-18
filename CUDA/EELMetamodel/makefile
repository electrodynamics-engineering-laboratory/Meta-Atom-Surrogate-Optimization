CC=g++
NVCC=nvcc
CXXFLAGS= -fopenmp -03 -Wextra -std=c++11
CUDAFLAGS= -std=c++11

all: metamodel

EELSurrogate.o: EELSurrogate.cuh EELSurrogate.cu
	nvcc $(CUDAFLAGS) -c EELSurrogate.cu

kernel.o: kernel.cu
	nvcc $(CUDAFLAGS) -c kernel.cu

metamodel: kernel.o EELSurrogate.o
	nvcc $(CUDAFLAGS) -o metamodel.out kernel.o EELSurrogate.o

clean:
	rm *.o *~
