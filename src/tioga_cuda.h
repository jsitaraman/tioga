#ifndef TIOGA_CUDA_H
#define TIOGA_CUDA_H

#include <cstdlib>
#include <cuda.h>

#define TIOGA_GPU_DEVICE __device__
#define TIOGA_GPU_GLOBAL __global__
#define TIOGA_GPU_HOST __host__
#define TIOGA_GPU_HOST_DEVICE __host__ __device__

namespace TIOGA {
namespace gpu {

using gpuError_t = cudaError_t;
constexpr gpuError_t gpuSuccess = cudaSuccess;

inline gpuError_t gpuGetLastError() { return cudaGetLastError(); }
inline const char* gpuGetErrorString(gpuError_t err) { return cudaGetErrorString(err); }

#define TIOGA_CUDA_CHECK_ERROR(call) {                                            \
        TIOGA::gpu::gpuError_t gpu_ierr = (call);                                 \
        if (TIOGA::gpu::gpuSuccess != gpu_ierr) {                                 \
            std::string err_str(std::string("TIOGA GPU error: ") + __FILE__       \
                               + ":" + std::to_string(__LINE__)                   \
                               + ": " + TIOGA::gpu::gpuGetErrorString(gpu_ierr)); \
            throw std::runtime_error(err_str);                                    \
        }}

template <typename T>
inline T* allocate_on_device(const size_t size)
{
    T* dptr = nullptr;
    TIOGA_CUDA_CHECK_ERROR(cudaMalloc((void**)(&dptr), size));
    return dptr;
}

template <typename T>
inline void copy_to_device(T* dptr, const T* hptr, const size_t size)
{
    TIOGA_CUDA_CHECK_ERROR(cudaMemcpy(dptr, hptr, size, cudaMemcpyHostToDevice));
}

template <typename T>
inline T* push_to_device(const T* hptr, const size_t size)
{
    T* dptr = allocate_on_device<T>(size);
    TIOGA_CUDA_CHECK_ERROR(cudaMemcpy(dptr, hptr, size, cudaMemcpyHostToDevice));
    return dptr;
}

template <typename T>
inline void pull_from_device(T* hptr, T* dptr, const size_t size)
{
    TIOGA_CUDA_CHECK_ERROR(cudaMemcpy(hptr, dptr, size, cudaMemcpyDeviceToHost));
}

template <typename T>
inline void deallocate_device(T** dptr)
{
    TIOGA_CUDA_CHECK_ERROR(cudaFree(static_cast<void*>(*dptr)));
    *dptr = nullptr;
}

} // namespace gpu
} // namespace TIOGA


#endif /* TIOGA_CUDA_H */
