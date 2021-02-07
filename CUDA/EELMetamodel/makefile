all: main

CC=g++
NVCC=nvcc
CXXFLAGS= -fopenmp -03 -Wextra -std=c++11
CUDAFLAGS= -std=c++11 -c

main: kernel.o surrogate.o
	nvcc $(CUDAFLAGS) kernel.o surrogate.o -o metamodel.out

surrogate.o: EELSurrogate.cuh EELSurrogate.cu
	nvcc $(CUDAFLAGS) EELSurrogate.cu -o surrogate.o

kernel.o: kernel.cu
	nvcc $(CUDAFLAGS) kernel.cu -o kernel.o

clean:
	rm *.o *~
