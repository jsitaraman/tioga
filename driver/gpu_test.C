#include <limits>
#include <vector>
#include <iostream>
#include "tioga_gpu.h"

void print_gpu_info()
{
    namespace gpu = TIOGA::gpu;
    std::cout << "BEGIN TEST print_gpu_info" << std::endl;
#if defined(TIOGA_HAS_GPU) && !defined(TIOGA_FAKE_GPU)
#if defined(CUDA_VERSION)
    std::cout << "CUDA configuration: "
              << "CUDA_VERSION: " << CUDA_VERSION
              << " " << CUDA_VERSION / 1000 << " "
              << (CUDA_VERSION % 1000) / 10 << std::endl;
#endif
    int ndevices;
    gpu::gpuDeviceProp_t dev;

    TIOGA_GPU_CALL_CHECK(GetDeviceCount(&ndevices));
    std::cout << "Num. devices = " << ndevices << std::endl;

    int rankDevice;
    TIOGA_GPU_CALL_CHECK(GetDevice(&rankDevice));
    TIOGA_GPU_CALL_CHECK(GetDeviceProperties(&dev, rankDevice));
    char busid[512];
    TIOGA_GPU_CALL_CHECK(DeviceGetPCIBusId(busid, 512, rankDevice));
    std::cout << "[" << rankDevice << "] "
              << dev.name << " CC: " << dev.major << "." << dev.minor
              << " ID: " << busid
              << " GM: "
              << (static_cast<double>(dev.totalGlobalMem) / (1 << 30)) << "GB"
              << " ShMem/Blk: " << (dev.sharedMemPerBlock / (1 << 10)) << "KB"
              << std::endl;

#else
    std::cout << "Not compiled with GPU support" << std::endl;
#endif
    std::cout << "END TEST print_gpu_info" << std::endl << std::endl;
}

void test_gpu_send_recv()
{
    namespace gpu = TIOGA::gpu;

    std::cout << "BEGIN TEST gpu_send_recv" << std::endl;
    std::vector<int> hvec(10, 10);

    int* dptr = gpu::push_to_device<int>(hvec.data(), hvec.size() * sizeof(int));

    if (dptr == nullptr) {
        std::cerr << "FAIL TEST gpu_send_recv" << std::endl;
        return;
    }

    gpu::memset_on_device(dptr, 0, hvec.size() * sizeof(int));

    gpu::pull_from_device(hvec.data(), dptr, hvec.size() * sizeof(int));

    int max_val = std::numeric_limits<int>::lowest();
    int min_val = std::numeric_limits<int>::max();
    for (int i=0; i < 10; ++i) {
        max_val = std::max(hvec[i], max_val);
        min_val = std::min(hvec[i], min_val);
    }

    std::cout << "Min = " << min_val << " Max = " << max_val << std::endl;

    TIOGA_FREE_DEVICE(dptr);

    if (dptr != nullptr)
        std::cout << "ERROR deallocating memory" << std::endl;
    std::cout << "END TEST gpu_send_recv" << std::endl << std::endl;
}

#ifdef TIOGA_HAS_GPU
TIOGA_GPU_GLOBAL void vec_add(double* out, double* a, double* b, int n)
{
    for (int i=0; i < n; ++i)
        out[i] = a[i] + b[i];
}

void test_cuda_vec_add()
{
    std::cout << "BEGIN TEST test_cuda_vec_add" << std::endl;
    namespace gpu = TIOGA::gpu;
    constexpr int N = 10;
    std::vector<double> avec(N, 10.0), bvec(N, 20.0), cvec(N, 0.0);

    double* da = gpu::push_to_device(avec.data(), avec.size() * sizeof(double));
    double* db = gpu::push_to_device(bvec.data(), bvec.size() * sizeof(double));
    double* dc = gpu::push_to_device(cvec.data(), cvec.size() * sizeof(double));
    TIOGA_GPU_LAUNCH_FUNC(vec_add, 1, 1, 0, 0, dc, da, db, N);

    gpu::pull_from_device(cvec.data(), dc, cvec.size() * sizeof(double));

    constexpr double csum_exact = 30.0 * N;
    double sum = 0.0;
    for (int i=0; i < N; ++i)
        sum += cvec[i];

    std::cout << "Expected = " << csum_exact << "; actual = " << sum << std::endl;

    TIOGA_FREE_DEVICE(da);
    TIOGA_FREE_DEVICE(db);
    TIOGA_FREE_DEVICE(dc);
    std::cout << "END TEST test_cuda_vec_add" << std::endl << std::endl;
}
#endif

int main()
{
    print_gpu_info();
    test_gpu_send_recv();

#ifdef TIOGA_HAS_GPU
    test_cuda_vec_add();
#endif

    return 0;
}
