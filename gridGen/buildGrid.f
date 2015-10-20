!
! Tioga is a library for overset grid assembly on parallel distributed systems
! Copyright (C) 2015 Jay Sitaraman
!
!
!     Create a unstructured grid partition
!     composed of prizm layers on a sphere and 
!     a Cartesian that overlaps it 
!
!     J. Sitaraman
!     04/18/2011
!
!     This code is not well commented and written
!     in partial f90 style.
!     Use with caution.
!     
!     begin main program
!
!
      program buildGrid

      use volumeData
      include 'mpif.h'
!
!     local variables
!
      character*128 fname
      integer :: myid,ierr,numprocs
      real*8  :: xlo(3),xhi(3),dx(3)

      namelist/inputs/nlayers,wallSpacing,outerSpacing,outerB,xlo,xhi,dx

      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,myid,ierr)
      call mpi_comm_size(mpi_comm_world,numprocs,ierr)

      open(unit=15,file='input',form='formatted')
      read(15,inputs)
      close(15)
      
      fname='data.tri'
      if (myid==0) write(6,*) 'Generating ',numprocs,' partitions' 
      !write(6,*) 'reading and constructing face info in partition :',myid
      if (myid < numprocs/2) then
       call getSurface(fname,myid,numprocs/2)
       call generatePrizms
      else
       nNodes=0
       nCells=0
       n6=0
       nwbc=0
       nobc=0
       call generateHex(myid-numprocs/2,xlo,xhi,dx,numprocs/2)
      endif
      !write(6,*) 'finished generating grid in partition :',myid
!
!      output volumes and be done!
!
      call mpi_barrier(mpi_comm_world,ierr)
      if (myid==0) write(6,*) 'Finished generating grid ..'
      if (myid==0) write(6,*) 'Writing tecplot compatible files..'
      call outputVolumes(myid)
      call mpi_barrier(mpi_comm_world,ierr)
      if (myid==0) write(6,*) 'All done'
      call mpi_finalize(ierr)
!
      end
!
!     
!     read face information
!
!
      subroutine getSurface(fname,myid,numprocs)

      use faceData

      implicit none
!
!     subroutine arguments
!
      character*128 fname
      integer :: myid,numprocs
!
!     local variables
!

      integer :: i,n,j,k
      integer :: nfacesTemp,nptsTemp,f(3)
      integer, allocatable :: itag(:)
      integer :: num,ib(3)
      real*8 :: xc(3),normTsq,vec1(3),vec2(3),vec3(3),arange(4)
      real*8 :: theta,phi,pi
!
!     begin
!
      pi=acos(-1.)
      open(unit=1,file=fname,form='formatted')

      read(1,*) nptsTemp,nfacesTemp

      allocate(itag(nptsTemp))
      allocate(xf(3,nptsTemp))
      allocate(normT(3,nptsTemp))

      do i=1,nptsTemp
         read(1,*) xf(1,i),xf(2,i),xf(3,i)
      enddo

      allocate(faceConn(3,nfacesTemp))

      num=myid
      do j=1,3
         ib(j)=mod(num,2)
         ib(j)=2*ib(j)-1
         num=num/2
      enddo
      !
      nfaces=0
      normT=0.
      !
      call getSpherePart(myid+1,numprocs,arange)
      !
      do i=1,nfacesTemp
         read(1,*) f(1),f(2),f(3)
         xc=0.
         do j=1,3
            do k=1,3
               xc(j)=xc(j)+xf(j,f(k))
            enddo
         enddo
         xc=xc/3
         theta=atan(xc(3)/sqrt(xc(2)**2+xc(1)**2))         
         phi=atan2(xc(2),xc(1))
         
         if (phi < 0) phi=phi+2*pi
         
         if ((theta-arange(1))*(theta-arange(2)) .le. 0. .and. 
     &        (phi-arange(3))*(phi-arange(4)).le.0) then
            nfaces=nfaces+1
            faceConn(1,nfaces)=f(1)
            faceConn(2,nfaces)=f(2)
            faceConn(3,nfaces)=f(3)
         endif

         vec1=xf(:,f(2))-xf(:,f(1))
         vec2=xf(:,f(3))-xf(:,f(1))

         vec3(1)=vec1(2)*vec2(3)-vec1(3)*vec2(2)
         vec3(2)=vec1(3)*vec2(1)-vec1(1)*vec2(3)
         vec3(3)=vec1(1)*vec2(2)-vec1(2)*vec2(1)
         
         do k=1,3
            normT(:,f(k))=normT(:,f(k))+vec3
         enddo

      enddo
      !
      itag=0
      do i=1,nfaces
         do j=1,3
            if (itag(faceConn(j,i))==0) then
               itag(faceConn(j,i))=1
            endif
         enddo
      enddo
      !
      npts=0
      !
      do i=1,nptsTemp
        if (itag(i) > 0) then
          npts=npts+1
          itag(i)=npts
        endif
      enddo
      !
      do i=1,nfaces
         do j=1,3
            faceConn(j,i)=itag(faceConn(j,i))
         enddo
      enddo
      !
      do i=1,nptsTemp
         if (itag(i) > 0) then
            xf(:,itag(i))=xf(:,i)
            normTsq=sqrt(normT(1,i)**2+
     &           normT(2,i)**2+
     &           normT(3,i)**2)

            normT(:,itag(i))=normT(:,i)/normTsq
         endif
      enddo

601   format('TITLE ="',a40,'"')
602   format('VARIABLES ="X", "Y", "Z"')
603   format('ZONE T="DATA",N=',i8,' E=',i8,' ET=QUADRILATERAL', 
     &  ' F=FEPOINT')

      ! 
      close(1)
      return
      end
!
!     generate Hex
!
      subroutine generateHex(myid,xl,xh,ddx,numprocs)
      use volumeData
      implicit none
      !
      ! subroutine arguments
      !
      integer :: myid,numprocs
      real*8  :: xl(3),xh(3),ddx(3)
      !
      ! local variables
      !
      integer :: i,j,k,l,nhpts,nhex,num,n
      integer ::ntot(3),nn(3)
      real*8 :: xlo(3),dx(3),dxx(3)
      real*8, allocatable :: x(:,:),tmpNode(:,:)
      integer, allocatable ::tmpCells(:,:),ndc8(:,:),tmpTag(:)
      integer :: nx(3)
      !
      call factorize(numprocs,nx)
      !
      do n=1,3
         ntot(n)=nint((xh(n)-xl(n))/ddx(n))
         num=mod(ntot(n)-1,nx(n))
         ntot(n)=max(ntot(n)-num,2)
         dx(n)=(xh(n)-xl(n))/(ntot(n)-1)
         dxx(n)=(ntot(n)-1)/nx(n)*dx(n)
         nn(n)=dxx(n)/dx(n)+1
      enddo
      !
      num=0
      do i=1,nx(1)
         do j=1,nx(2)
            do k=1,nx(3)
               num=num+1
               if (myid+1 == num) then
                  xlo(1)=xl(1)+(i-1)*dxx(1)
                  xlo(2)=xl(2)+(j-1)*dxx(2)
                  xlo(3)=xl(3)+(k-1)*dxx(3)
               endif
            enddo
         enddo
      enddo
      !
      !write(6,*) 'myid,nx=',myid,nx
      !write(6,*) 'hex grid in proc :',myid,'nn=',nn,'dx=',dx
      !
      ! allocate memory for grid and connectivity
      !
      allocate(x(3,nn(1)*nn(2)*nn(3)))
      allocate(ndc8(8,(nn(1)-1)*(nn(2)-1)*(nn(3)-1)))
      !
      ! create coordinate data
      !
      nhpts=0
      do l=1,nn(3)
         do k=1,nn(2)
            do j=1,nn(1)
               nhpts=nhpts+1
               x(1,nhpts)=(j-1)*dx(1)+xlo(1)
               x(2,nhpts)=(k-1)*dx(2)+xlo(2)
               x(3,nhpts)=(l-1)*dx(3)+xlo(3)
            enddo
         enddo
      enddo
      !
      ! create hex connectivity 
      !
      nhex=0
      do l=1,nn(3)-1
         do k=1,nn(2)-1
            do j=1,nn(1)-1
               nhex=nhex+1
               ndc8(1,nhex)=(l-1)*nn(2)*nn(1)+(k-1)*nn(1)+j
               ndc8(2,nhex)=ndc8(1,nhex)+1
               ndc8(3,nhex)=ndc8(2,nhex)+nn(1)
               ndc8(4,nhex)=ndc8(1,nhex)+nn(1)
               ndc8(5:8,nhex)=ndc8(1:4,nhex)+nn(1)*nn(2)
            enddo
         enddo
      enddo
      !
      ! now append hexes to the original list of prizms
      !
      n6=nCells
      n8=nhex
      if (nNodes > 0) then
       allocate(tmpNode(3,nNodes),tmpCells(8,nCells),tmpTag(nNodes))
       !
       do i=1,nCells
          tmpCells(:,i)=volConn(:,i)
       enddo
       do i=1,nNodes
          tmpNode(:,i)=xv(:,i)
          tmpTag(i)=bodyTag(i)
       enddo
       !
       deallocate(xv,volConn,bodyTag)
      endif
      !
      nCells=nCells+nhex
      nNodes=nNodes+nhpts
      allocate(volConn(8,nCells),xv(3,nNodes))
      allocate(bodyTag(nNodes))
      !
      if (nCells-nhex > 0) then
       do i=1,nCells-nhex
         volConn(:,i)=tmpCells(:,i)
       enddo
      endif
      !
      do i=1,nhex
         volConn(:,i+nCells-nhex)=ndc8(:,i)+nNodes-nhpts
      enddo
      !
      if (nNodes-nhpts > 0) then
       do i=1,nNodes-nhpts
         xv(:,i)=tmpNode(:,i)
         bodyTag(i)=tmpTag(i)
       enddo
      endif
      !
      do i=1,nhpts
         xv(:,i+nNodes-nhpts)=x(:,i)
         bodyTag(i+nNodes-nhpts)=2
      enddo
      !
      if (allocated(tmpNode)) deallocate(tmpNode)
      if (allocated(tmpTag)) deallocate(tmpTag)
      if (allocated(tmpCells)) deallocate(tmpCells)
      !
      deallocate(x,ndc8)
      !
      return
      end          
!     
!     
!     generate  Prizms
!
!
      subroutine generatePrizms

      use volumeData
      use faceData

      implicit none
!
!     local variables
!
      integer :: l,f,n,a,b,c,a1,b1,c1,n1,i
      integer :: lsum,lcells,na,nb
      real *8 :: cntrlpt(2),arclen(2),ds(2)
      real *8 :: vec1(3),vec2(3),vec3(3)
      real *8 :: normTsq,ss
      integer, allocatable :: inorm(:)
      real *8, allocatable :: strandSpacing(:)
      real *8, allocatable :: faceConnT(:,:)
!
!     begin execution
!
!
!     allocate global arrays
!
      nNodes=npts*(nlayers+1)
      nCells=nfaces*nlayers
      n6=nCells
      allocate(xv(3,nnodes))
      allocate(volConn(8,nCells))
!
!     allocate local arrays
!     
      allocate(faceConnT(3,nfaces))
      allocate(strandSpacing(nlayers+1))
!
!     create hyperbolic distribution in strand
!
      cntrlpt(1)=1
      ds(1)=wallSpacing
      arclen(1)=0.

      cntrlpt(2)=nlayers+1
      ds(2)=outerSpacing
      arclen(2)=outerB

      !write(6,*) ds
      call distrib(2,cntrlpt,arclen,ds,strandSpacing)
      !do l=1,nlayers+1
      !write(6,*) strandSpacing(l)
      !enddo
!
!     initalize the marching boundary
!      
      faceConnT(:,1:nfaces)=faceConn(:,1:nfaces)
      xv(:,1:npts)=xf(:,1:npts)
      lsum=0
      lcells=0
!
!     begin generating volumes
!    
      do l=1,nlayers
         na=lsum+1
         nb=lsum+npts
         do n=na,nb
            n1=n-lsum
            ss=strandSpacing(l+1)-strandSpacing(l)
            xv(:,npts+n)=xv(:,n)+normT(:,n1)*ss
         enddo
         do f=1,nfaces
            volConn(1:3,lcells+f)=faceConnT(:,f)
            volConn(4,lcells+f)=volConn(3,lcells+f)
            volConn(5:7,lcells+f)=faceConnT(:,f)+npts
            volConn(8,lcells+f)=volConn(7,lcells+f)
         enddo
         lsum=lsum+npts
         faceConnT=faceConnT+npts
         lcells=lcells+nfaces
      enddo

      nwbc=npts
      nobc=npts
      !
      ! set wall and overset bc indices
      ! because of simple geometry here
      ! the first npts nodes are wall nodes
      ! last npts nodes are overset bc nodes
      !
      allocate(wbcnode(nwbc))
      allocate(obcnode(nwbc))
      !
      do i=1,nwbc
         wbcnode(i)=i
      enddo
      !
      do i=1,nobc
         obcnode(i)=nNodes-nobc+i
      enddo
      !
      allocate(bodyTag(nNodes))
      do i=1,nNodes
         bodyTag(i)=1
      enddo
      !
      deallocate(strandSpacing)
      !
      return
      end
! 
!        
!     write tecplot compatible output file 
!
      subroutine outputVolumes(myid)
!
      use volumeData 
      implicit none
!
      integer myid
      integer i,n
      character*128 :: fname,integer_string

      write(integer_string,"(I7)") myid
      
      fname='cell'//trim(adjustl(integer_string))//'.plt'
      !write(6,*) 'writing partition :',myid,'to ',trim(adjustl(fname))
      open(unit=10,file=fname,form='formatted')
      write(10,601) n6,n8,nNodes,nCells,nwbc,nobc
      write(10,*) 'TITLE =" Unstructured grid"   '
      write(10,*) 'VARIABLES ="X", "Y", "Z", "bodyTag"'
      write(10,606) nNodes,nCells
      
      do i=1,nnodes
       write(10,605) xv(1,i),xv(2,i),xv(3,i),bodyTag(i)
      enddo

      do i=1,ncells
       write(10,1004) volConn(1,i),volConn(2,i),volConn(3,i),
     &              volConn(4,i),volConn(5,i),volConn(6,i),
     &               volConn(7,i),volConn(8,i)
      enddo
      do i=1,nwbc
         write(10,*) wbcnode(i)
      enddo
      do i=1,nobc
         write(10,*) obcnode(i)
      enddo

      close(10)
  601 format('# ',6(1X,I10))
  606 format('ZONE T="VOL_MIXED",N=',i8,' E=',i8,' ET=BRICK,',
     .' F=FEPOINT')
  605 format((3(E15.8,1X),I7))
           
1004  format(8(1X,I10))
      return   
      end

