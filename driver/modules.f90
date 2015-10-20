!
! Tioga is a library for overset grid assembly on parallel distributed systems
! Copyright (C) 2015 Jay Sitaraman
!
!
module gridtype
 type grid
    integer :: nv,n4,n5,n6,n8,nwbc,nobc,nvar,ncells,nghost,nmax,ndof
    real*8, allocatable :: x(:),q(:),dq(:),xcentroid(:)
    integer, allocatable :: bodytag(:)
    integer, allocatable :: iblank(:)
    integer, allocatable :: ghostData(:,:)
    integer, allocatable :: ndc4(:,:),ndc5(:,:),ndc6(:,:),ndc8(:,:)
    integer, allocatable :: wbcnode(:),obcnode(:)
    real*8, allocatable :: scal(:)
 end type grid
end module gridtype
