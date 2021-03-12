#ifndef TIOGA_NOGPU_H
#define TIOGA_NOGPU_H

#include <cstdlib>
#include <cstring>

#define TIOGA_GPU_DEVICE
#define TIOGA_GPU_GLOBAL
#define TIOGA_GPU_HOST
#define TIOGA_GPU_HOST_DEVICE

namespace TIOGA {
namespace gpu {

using gpuError_t = int;
constexpr gpuError_t gpuSuccess = 0;

inline gpuError_t gpuGetLastError() { return gpuSuccess; }
inline const char* gpuGetErrorString(gpuError_t err) { return "Success"; }

#define TIOGA_GPU_LAUNCH_FUNC(func, blocks, threads, sharedmem, stream, ...) \
  func(__VA_ARGS__);

template <typename T>
inline T* allocate_on_device(const size_t size)
{
    return static_cast<T*>(std::malloc(size));
}

template <typename T>
inline void copy_to_device(T* dptr, const T* hptr, const size_t size)
{
    std::memcpy(dptr, hptr, size);
}

template <typename T>
inline T* push_to_device(const T* hptr, const size_t size)
{
    T* dptr = allocate_on_device<T>(size);
    std::memcpy(dptr, hptr, size);
    return dptr;
}

template <typename T>
inline void pull_from_device(T* hptr, T* dptr, const size_t size)
{
    std::memcpy(hptr, dptr, size);
}

template <typename T>
inline void deallocate_device(T** dptr)
{
    std::free(*dptr);
    *dptr = nullptr;
}

template <typename T>
inline void memset_on_device(T* dptr, T val, const size_t sz)
{
  std::memset(dptr, val, sz);
}

inline void synchronize(void)
{}
}
}


#endif /* TIOGA_NOGPU_H */
