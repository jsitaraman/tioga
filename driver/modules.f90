!
! Tioga is a library for overset grid assembly on parallel distributed systems
! Copyright (C) 2015 Jay Sitaraman
!
!
module gridtype
 type grid
    integer :: nv,n4,n5,n6,n8,nwbc,nobc,nvar,ncells,nghost,nmax,ndof
    integer :: mdim(3)
    real*8  :: dx(3),dx0(3)
    real*8,  pointer :: x(:),q(:),s(:),dq(:),xcentroid(:),q0(:)
    integer, allocatable :: bodytag(:)
    integer, allocatable :: iblank(:)
    integer, allocatable :: ghostData(:,:)
    integer, allocatable :: ndc4(:,:),ndc5(:,:),ndc6(:,:),ndc8(:,:)
    integer, allocatable :: wbcnode(:),obcnode(:)
    real*8, allocatable :: scal(:)
    real*8, allocatable :: nres(:),cres(:)
 end type grid
end module gridtype
