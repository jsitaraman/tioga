!
! Tioga is a library for overset grid assembly on parallel distributed systems
! Copyright (C) 2015 Jay Sitaraman
!
!> Poisson MMS test code  for Topology Independent Overset Grid Assembler (TIOGA)
!> 
!> Jay Sitaraman
!>
!> 03/05/2014
!>
program poissonMMS
  !
  use gridtype
  !
  implicit none
  !
  include 'mpif.h'
  !
  integer, parameter :: mgrids=4
  type(grid), target :: gr(mgrids)
  type(grid), pointer :: g
  !
  integer :: myid,numprocs,ierr
  integer :: ib,iter,ngrids,irefine,nref,ndof,ntypes,nv2,i,jmax,jmax2
  integer :: mexclude,nfringe
  integer :: mdim(3)
  !
  logical :: testSolver,saveout
  !
  real*8  :: BIGVALUE,SMALLVALUE,x,z,overlap,TOL,res,err,rfac
  real*8 :: slope(mgrids),dx,bbox_inner(2),bbox_outer(2)
  real*8 :: xhi(3),xlo(3),e0(mgrids),e(mgrids)
  !
  character*32 :: refstr,rfacstr,overlapstr,saveoutstr,testSolverstr
  !
  ! initialize mpi
  !
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,myid,ierr)
  call mpi_comm_size(mpi_comm_world,numprocs,ierr)
  !
  ! initialize grids
  !
  ngrids=2            ! number of grids (only works for 2 in this example)
  !
  testSolver=.false.  ! .true. does Gauss-Seidel iterations, otherwise just tests a laplacian against
                      ! exact
  saveout = .false.   ! save output grids
  !
  nref=5                ! number of refinement cycles
  rfac=2d0              ! factor for refinement
  overlap=0d0           ! physical overlap distance between fringes of coarse and fine grids
  !
  ! just a basic argument parser
  ! need to improve this
  !
  if (iargc() > 0) then
    if (iargc()>=1) then
        call getarg(1,refstr)
        read(refstr,*) nref
     endif
     if (iargc() >=2) then
        call getarg(2,rfacstr)
        read(rfacstr,*) rfac
     endif
     if (iargc() >=3) then
        call getarg(3,overlapstr)
        read(overlapstr,*) overlap
     endif
     if (iargc() >=4) then
        call getarg(4,saveoutstr)
        read(saveoutstr,*) saveout
     endif
     if (iargc() >=5) then
        call getarg(5,testSolverstr)
        read(testSolverstr,*) testSolver
     endif
  else
     write(6,*) 'poisson_mms <nref> <rfac> <overlap> <saveout> <testSolver>'
     write(6,*) '------- using default parameters ----------------'
  endif
  !
  write(6,*) '----------------------------------'
  write(6,"(A20,1x,I10)") 'nref:',nref
  write(6,"(A20,1x,F10.4)") 'rfac:',rfac
  write(6,"(A20,1x,F10.4)") 'overlap:',overlap
  write(6,"(A20,1x,L10)") 'savout:',saveout
  write(6,"(A20,1x,L10)") 'testSolver:',testSolver
  write(6,*) '----------------------------------'
  !
  ! bounding boxes to be differentiate regions for creating node and
  ! cell resolution fields
  !
  bbox_outer=[0.25d0,0.75d0]
  bbox_inner=[0.25d0-overlap,0.75d0+overlap]
  !
  ! resolutions for different regions
  !
  BIGVALUE=1d15 
  SMALLVALUE=2d0 
  TOL=1e-5
  !
  ! mexclude and nfringe for
  ! this problem
  !
  mexclude=1
  nfringe=1
  !
  open(unit=101,file='error.dat',status='unknown')
  !
  do irefine=1,nref
     !
     ! intialize tioga
     !
     call tioga_init_f90(mpi_comm_world)
     call tioga_setnfringe(nfringe)
     call tioga_setmexclude(mexclude)
     !
     ! create the grids
     !
     do ib =1, ngrids
        if (ib==1) then
           xlo=[-0.5d0,   0d0, -0.5d0]     ! hard coded corners of coarse grid
           xhi=[ 1.5d0, 1d0/16,  1.5d0]
           jmax=rfac**(irefine-1)*32+1
           mdim=[jmax,3,jmax]
        else if (ib==2) then
           xlo=[0d0,      0d0, 0d0]      ! hard coded corners of fine grid
           xhi=[1d0, 1d0/16d0, 1d0]
           jmax2=rfac**(irefine-1)*28+1
           mdim=[jmax2,3,jmax2]
        endif
        if (irefine > 1) call deallocategrid(gr(ib))
        call initgrid(gr(ib),xlo,xhi,mdim)
     enddo
     !
     ! register grid data
     !
     ntypes=1
     nv2=8
     do ib=1,ngrids
        g=>gr(ib)
        g%nobc=0  !-> add this to prevent mandatory receptor rules being applied
        call tioga_registergrid_data_mb(ib,ib,g%nv,g%x,g%iblank,&
             g%nwbc,g%nobc,g%wbcnode,g%obcnode,&
             ntypes,nv2,g%n8,g%ndc8)
     enddo
     !
     ! set the resolutions
     !     
     do ib=1,ngrids
        g=>gr(ib)
        allocate(g%nres(g%nv))
        allocate(g%cres(g%n8))
        do i=1,g%nv
           x=g%x(3*i-2);z=g%x(3*i)
           if (ib==1) then
              if (x > bbox_outer(1)-TOL .and. x < bbox_outer(2)+TOL .and. &
                  z > bbox_outer(1)-TOL .and. z < bbox_outer(2)+TOL) then
                 g%nres(i)=BIGVALUE
              else
                 g%nres(i)=SMALLVALUE
              endif
           elseif (ib==2) then
              if (x > bbox_inner(1)+TOL .and. x < bbox_inner(2)-TOL .and. &
                  z > bbox_inner(1)+TOL .and. z < bbox_inner(2)-TOL) then
                 g%nres(i)=SMALLVALUE
              else
                 g%nres(i)=BIGVALUE
              endif
           endif
        enddo
        !
        do i=1,g%n8
           x=g%xcentroid(3*i-2)
           z=g%xcentroid(3*i)
           if (x > bbox_inner(1)+TOL .and. x < bbox_inner(2)-TOL .and. &
                z > bbox_inner(1)+TOL .and. z < bbox_inner(2)-TOL) then
              if (ib==1) g%cres(i)=BIGVALUE
              if (ib==2) g%cres(i)=SMALLVALUE
           else
              if (ib==1) g%cres(i)=SMALLVALUE
              if (ib==2) g%cres(i)=BIGVALUE
           endif
        enddo
        call tioga_setresolutions_multi(ib,g%nres,g%cres)
     end do
     !
     call tioga_preprocess_grids
     call tioga_performconnectivity
     !
     ! solve poisson MMS
     !
     if (testSolver) then
        res=1d0
        iter=0
        do while(res > 1e-14)
           iter=iter+1
           do ib=1,ngrids
              call tioga_registersolution(ib,gr(ib)%q)
           enddo
           do ib=1,ngrids
              call tioga_dataupdate_mb(gr(ib)%nvar,'row')
           enddo
           !
           ! gauss-seidel iteration
           !
           do ib=1,ngrids
              call iterate(gr(ib),res,1)
           enddo
           do ib=1,ngrids
              call tioga_dataupdate_mb(gr(ib)%nvar,'row')
           enddo
           if (mod(iter,100) ==0) write(6,*) iter,res,ib
        end do
        !
        ! compute algebraic error
        !
        do ib=1,ngrids
           call compute_alg_err(gr(ib),e(ib))
           if (irefine > 1) then
              slope(ib)=(log(e0(ib))-log(e(ib)))/(log(gr(ib)%dx0(1))-log(gr(ib)%dx(1)))
           endif
           gr(ib)%dx0=gr(ib)%dx
        enddo
     else
        !
        ! check the laplacian operator against exact
        !
        do ib=1,ngrids
           call tioga_registersolution(ib,gr(ib)%q)
        enddo
        do ib=1,ngrids
           call tioga_dataupdate_mb(gr(ib)%nvar,'row')
        enddo
        do ib=1,ngrids
           call poisson_operator(gr(ib),e(ib))
           if (irefine > 1) then
              slope(ib)=(log(e0(ib))-log(e(ib)))/(log(gr(ib)%dx0(1))-log(gr(ib)%dx(1)))
           endif
           gr(ib)%dx0=gr(ib)%dx
        enddo
     endif
     !
     if (irefine > 1 ) then
        write(6,1000) jmax,jmax2,e(1:ngrids),slope(1:ngrids)
        write(101,1000) jmax,jmax2,e(1:ngrids),slope(1:ngrids)
     else
        write(6,1000) jmax,jmax2,e(1:ngrids)
     endif
     e0=e
     if (saveout) call tioga_writeoutputfiles(g%nvar,'row') !< write output files, if need be
     call tioga_delete
  enddo
  !
  close(101)
  !
1000 format((I5,1x,I5,1x,12(2x,F10.4)))
  !
  call mpi_finalize(ierr)
  !
end program poissonMMS
