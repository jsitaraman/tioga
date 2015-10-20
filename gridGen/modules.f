!
! Tioga is a library for overset grid assembly on parallel distributed systems
! Copyright (C) 2015 Jay Sitaraman
!
!
!     surface data
!
      module faceData

      integer :: npts,nfaces                !< number of nodes and faces
      real, allocatable :: xf(:,:)          !< list of nodes
      integer, allocatable :: faceConn(:,:) !< face connectivity
      real, allocatable :: normT(:,:)       !< vertex normals
 
      end module faceData
!
!     volume data
!
      module volumeData

      integer :: nCells,nNodes,nwbc,nobc               !< total number of cells, nodes
                                                       !< wall boundary nodes, overset
                                                       !< bc nodes
      integer :: n6,n8                                 !< number of prizm and hex


      real, allocatable :: xv(:,:)                     !< list of vertices
      integer, allocatable :: volConn(:,:)             !< prizm connectivity
      integer, allocatable :: hexConn(:,:)             !< hex connectivity
      integer, allocatable :: wbcnode(:)               !< list of wall nodes
      integer, allocatable :: obcnode(:)               !< list of obc nodes
      integer, allocatable :: bodyTag(:)               !< bodytag
      real *8 :: wallSpacing,outerSpacing,outerB,vfac  !< parameters for builing 
                                                       !< prizm layers
       integer :: nlayers
       !data nlayers /60/
       !data vfac /1.1/
       !data wallSpacing,outerSpacing,outerB /0.02,0.1,10.0/

      end module volumeData
