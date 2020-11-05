module kokkos
  use, intrinsic :: iso_c_binding

  implicit none

  public

    interface
       subroutine f_kokkos_initialize() &
           bind(c, name='c_kokkos_initialize')
         use, intrinsic :: iso_c_binding
         implicit none
       end subroutine f_kokkos_initialize
    end interface

    interface
       subroutine f_kokkos_finalize() &
           bind(c, name='c_kokkos_finalize')
         use, intrinsic :: iso_c_binding
         implicit none
       end subroutine f_kokkos_finalize
    end interface

    contains

      subroutine kokkos_initialize()
        use, intrinsic :: iso_c_binding
        implicit none
        call f_kokkos_initialize()
      end subroutine kokkos_initialize

      subroutine kokkos_finalize()
        use, intrinsic :: iso_c_binding
        implicit none
        call f_kokkos_finalize()
      end subroutine kokkos_finalize

end module kokkos
