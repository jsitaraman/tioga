#ifndef TIOGA_GPU_H
#define TIOGA_GPU_H

#if defined(TIOGA_HAS_CUDA)
#include "tioga_cuda.h"
#elif defined(TIOGA_HAS_HIP)
#include "tioga_hip.h"
#else
#include "tioga_nogpu.h"
#endif

namespace TIOGA {
namespace gpu {

#if defined(TIOGA_HAS_GPU) && !defined(TIOGA_FAKE_GPU)
#define TIOGA_GPU_CHECK_ERROR(call) {                                             \
        TIOGA::gpu::gpuError_t gpu_ierr = (call);                                 \
        if (TIOGA::gpu::gpuSuccess != gpu_ierr) {                                 \
            std::string errStr(std::string("TIOGA GPU error: ") + __FILE__        \
                               + ":" + std::to_string(__LINE__)                   \
                               + ": " + TIOGA::gpu::gpuGetErrorString(gpu_ierr)); \
            throw std::runtime_error(errStr);                                     \
        }}
#else
#define TIOGA_GPU_CHECK_ERROR(call)  (call)
#endif

#define TIOGA_FREE_DEVICE(dptr) TIOGA::gpu::deallocate_device(&dptr)

}
}


#endif /* TIOGA_GPU_H */
