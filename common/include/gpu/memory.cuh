#pragma once
#include <cuda_runtime.h>
#include <stdio.h>

template <class T>
void MallocDeviceMemory(T **ptr, size_t n)
{
    cudaError_t err = cudaMalloc((void**)ptr, n * sizeof(T));
    if (err != cudaSuccess) {
        fprintf(stderr, "Error allocating device memory: %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

template <class T>
void FreeDeviceMemory(T *ptr)
{
    cudaError_t err = cudaFree(ptr);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error freeing device memory: %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

template <class T>
void CopyToDevice(T *dst, const T *src, size_t n)
{
    cudaError_t err = cudaMemcpy(dst, src, n * sizeof(T), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "Error copying to device memory: %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}