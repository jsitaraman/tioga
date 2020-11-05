#include <Kokkos_Core.hpp>

extern "C"
{
  void c_kokkos_initialize() { Kokkos::initialize(); }
  void c_kokkos_finalize() { Kokkos::finalize(); }
}
