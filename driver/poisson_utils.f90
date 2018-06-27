!
! Tioga is a library for overset grid assembly on parallel distributed systems
! Copyright (C) 2015 Jay Sitaraman
!
!=====================================================================
!
! initialize a cartesian grid
!
!=====================================================================
subroutine initgrid(g,xlo,xhi,mdim)
!
use gridtype
implicit none
!
type(grid), intent(inout) :: g
real*8, intent(in) :: xlo(3),xhi(3)
integer, intent(in) :: mdim(3)
!
integer :: i,j,k,l,n,m
real*8 :: pi,xx,yy,zz
!
! local variables
!
pi=acos(-1d0)
g%n4=0
g%n5=0
g%n6=0
g%n8=1
g%nv=1
do n=1,3
  g%mdim(n)=mdim(n)
  g%dx(n)=(xhi(n)-xlo(n))/(mdim(n)-1)
  g%n8=g%n8*(mdim(n)-1)
  g%nv=g%nv*mdim(n)
end do
!
g%nvar=1
g%ncells=g%n8
g%nghost=0
g%ndof=g%ncells+g%nghost
!
allocate(g%x(3*g%nv),g%s(g%nvar*g%nv),g%bodytag(g%nv),g%iblank(g%nv))
allocate(g%scal(g%nvar))
allocate(g%xcentroid(3*g%ndof))
allocate(g%q(g%nvar*g%nv))
allocate(g%q0(g%nvar*g%nv))
allocate(g%dq(3*g%nvar*g%nv))
allocate(g%ndc8(8,g%n8))
!
g%iblank=1
!
! read coordinates and body tag 
! also create some conservative variables
!
n=0
m=0
do l=1,mdim(3)
 zz=xlo(3)+g%dx(3)*(l-1)
 do k=1,mdim(2)
   yy=xlo(2)+g%dx(2)*(k-1)
   do j=1,mdim(1)
     xx=xlo(1)+g%dx(1)*(j-1)
     n=n+1
     m=m+1
     g%x(n)=xx
     g%q(m)=cos(2*pi*xx)*cos(2*pi*zz) !cos(2*pi*yy) !*cos(2*pi*zz)
     g%q0(m)=g%q(m)
     g%s(m)=-2*4d0*(pi**2d0)*g%q(m)
     n=n+1
     g%x(n)=yy
     n=n+1
     g%x(n)=zz
    enddo
 enddo
enddo
n=0
do l=1,mdim(3)-1
 do k=1,mdim(2)-1
  do j=1,mdim(1)-1
    n=n+1
    g%ndc8(1,n)=(l-1)*mdim(2)*mdim(1)+(k-1)*mdim(1)+j
    g%ndc8(2,n)=g%ndc8(1,n)+1
    g%ndc8(3,n)=g%ndc8(2,n)+mdim(1)
    g%ndc8(4,n)=g%ndc8(1,n)+mdim(1)
    g%ndc8(5:8,n)=g%ndc8(1:4,n)+mdim(2)*mdim(1)
  enddo
 enddo
enddo
!
m=0
do i=1,g%n8
    m=m+1
    g%xcentroid(3*m-2:3*m)=0.
    do j=1,8
       g%xcentroid(3*m-2)=g%xcentroid(3*m-2)+g%x(3*g%ndc8(j,i)-2)/8d0
       g%xcentroid(3*m-1)=g%xcentroid(3*m-1)+g%x(3*g%ndc8(j,i)-1)/8d0
       g%xcentroid(3*m)=g%xcentroid(3*m)+g%x(3*g%ndc8(j,i))/8d0
    enddo
enddo
!
! read wall and overset boundary nodes
!
g%nwbc=0
g%nobc=2*g%mdim(1)*g%mdim(2)+2*(g%mdim(2)*(g%mdim(3)-2)) !+2*(g%mdim(3)-2)*(g%mdim(1)-2)
allocate(g%wbcnode(g%nwbc),g%obcnode(g%nobc))
!
n=0
l=1
do k=1,mdim(2)
 do j=1,mdim(1)
  n=n+1
  g%obcnode(n)=(l-1)*mdim(2)*mdim(1)+(k-1)*mdim(1)+j
 enddo
enddo
l=mdim(3)
do k=1,mdim(2)
 do j=1,mdim(1)
  n=n+1
  g%obcnode(n)=(l-1)*mdim(2)*mdim(1)+(k-1)*mdim(1)+j
 enddo
enddo
j=1
do l=2,mdim(3)-1
 do k=1,mdim(2)
  n=n+1
  g%obcnode(n)=(l-1)*mdim(2)*mdim(1)+(k-1)*mdim(1)+j
 enddo
enddo
j=mdim(1)
do l=2,mdim(3)-1
 do k=1,mdim(2)
  n=n+1
  g%obcnode(n)=(l-1)*mdim(2)*mdim(1)+(k-1)*mdim(1)+j
 enddo
enddo
!k=1
!do j=2,mdim(1)-1
! do l=2,mdim(3)-1
!  n=n+1
!  g%obcnode(n)=(l-1)*mdim(2)*mdim(1)+(k-1)*mdim(1)+j
! enddo
!enddo
!k=mdim(2)
!do j=2,mdim(1)-1
! do l=2,mdim(3)-1
!  n=n+1
!  g%obcnode(n)=(l-1)*mdim(2)*mdim(1)+(k-1)*mdim(1)+j
! enddo
!enddo
!
return
!
end subroutine initGrid
!
! basic poisson operator
! 
subroutine poisson_operator(g,res)
!
use gridtype
implicit none
type(grid) :: g
real*8, intent(inout) :: res
!
integer :: i,j,k,l,idx
integer :: mdim(3)
real*8  :: qnew,diag,dx(3)
real*8, pointer :: q(:),s(:)
real*8, allocatable :: ee(:)
real*8 :: emax
integer :: indx(3)
!
mdim=g%mdim
q=>g%q
s=>g%s
dx=g%dx
!diag=2d0/dx(1)/dx(1)+2d0/dx(2)/dx(2)+2d0/dx(3)/dx(3)
diag=2d0/dx(1)/dx(1)+2d0/dx(3)/dx(3)
res=0d0
allocate(ee(mdim(1)*mdim(2)*mdim(3)))
!
ee=0d0
do l=2,mdim(3)-1
  do k=2,mdim(2)-1
    do j=2,mdim(1)-1
      idx=(l-1)*mdim(2)*mdim(1)+(k-1)*mdim(1)+j
      if (g%iblank(idx) .ne. 1) cycle
      qnew=  (q(idx+1)              + q(idx-1              ))/dx(1)/dx(1)+ &
             !(q(idx+mdim(1))        + q(idx-mdim(1)        ))/dx(2)/dx(2) !+ &
             (q(idx+mdim(1)*mdim(2))+ q(idx-mdim(1)*mdim(2)))/dx(3)/dx(3)
      res=res+(qnew-diag*q(idx)-s(idx))**2
      ee(idx)=abs(qnew-diag*q(idx)-s(idx))
      if (ee(idx) > emax) then
        emax=ee(idx)
        indx=[j,k,l]
      endif
    enddo
   enddo
enddo
do l=1,mdim(3)
  do k=1,mdim(2)
     do j=1,mdim(1)
      idx=(l-1)*mdim(2)*mdim(1)+(k-1)*mdim(1)+j
      q(idx)=ee(idx)
     enddo
  enddo
enddo
!
deallocate(ee)
res=sqrt(res/(mdim(3)-2)/(mdim(2)-2)/(mdim(1)-2))
l=indx(3)
j=indx(1)
idx=(l-1)*mdim(2)*mdim(1)+(k-1)*mdim(1)+j
!
return
end subroutine poisson_operator
!
! basic Gauss-Seidel Poisson solver
! 
subroutine iterate(g,res,niter)
!
use gridtype
implicit none
type(grid) :: g
real*8, intent(inout) :: res
integer,intent(in) :: niter
!
integer :: i,j,k,l,idx
integer :: mdim(3)
real*8  :: qnew,diag,dx(3)
real*8, pointer :: q(:),s(:)
!
mdim=g%mdim
q=>g%q
s=>g%s
dx=g%dx
res=0d0
!diag=2d0/dx(1)/dx(1)+2d0/dx(2)/dx(2)+2d0/dx(3)/dx(3)
diag=2d0/dx(1)/dx(1)+2d0/dx(3)/dx(3)
!
do i=1,niter
 do l=2,mdim(3)-1
  do k=2,mdim(2)-1
    do j=2,mdim(1)-1
      idx=(l-1)*mdim(2)*mdim(1)+(k-1)*mdim(1)+j
      if (g%iblank(idx) .ne. 1) cycle
      qnew=  (q(idx+1)              + q(idx-1              ))/dx(1)/dx(1)+ &
             !(q(idx+mdim(1))        + q(idx-mdim(1)        ))/dx(2)/dx(2)+ &
             (q(idx+mdim(1)*mdim(2))+ q(idx-mdim(1)*mdim(2)))/dx(3)/dx(3)
      qnew=(qnew-s(idx))/diag
      res=res+(qnew-q(idx))**2
      q(idx)=qnew
    enddo
   enddo
  enddo
enddo
res=sqrt(res/(mdim(3)-1)/(mdim(2)-1)/(mdim(1)-1))
!
return
end subroutine iterate

subroutine compute_alg_err(g,err)
!
use gridtype
implicit none
type(grid) :: g
real*8, intent(inout) :: err
!
integer :: i,j,k,l,mdim(3),idx,ecount
real*8 ::ee
!
err=0d0
ecount=0
mdim=g%mdim
do l=1,mdim(3)
  do k=1,mdim(2)
    do j=1,mdim(1)
       idx=(l-1)*mdim(2)*mdim(1)+(k-1)*mdim(1)+j
       if (g%iblank(idx).ne.1) cycle
       ee=abs(g%q0(idx)-g%q(idx))
       err=err+(g%q0(idx)-g%q(idx))**2
       g%q(idx)=ee
       ecount=ecount+1
    enddo
  enddo
enddo 
!err=sqrt(err/g%mdim(3)/g%mdim(2)/g%mdim(1))
err=sqrt(err/ecount)
end subroutine compute_alg_err

subroutine deallocategrid(g)
 use gridtype
 implicit none
 type(grid) :: g
 if (associated(g%x)) deallocate(g%x)
 if (associated(g%q)) deallocate(g%q)
 if (associated(g%s)) deallocate(g%s)
 if (associated(g%dq)) deallocate(g%dq)
 if (associated(g%xcentroid)) deallocate(g%xcentroid)
 if (associated(g%q0)) deallocate(g%q0)
 if (allocated(g%bodytag)) deallocate(g%bodytag)
 if (allocated(g%iblank)) deallocate(g%iblank)
 if (allocated(g%ghostData)) deallocate(g%ghostData)
 if (allocated(g%ndc4)) deallocate(g%ndc4)
 if (allocated(g%ndc5)) deallocate(g%ndc5)
 if (allocated(g%ndc6)) deallocate(g%ndc6)
 if (allocated(g%ndc8)) deallocate(g%ndc8)
 if (allocated(g%wbcnode)) deallocate(g%wbcnode)
 if (allocated(g%obcnode)) deallocate(g%obcnode)
 if (allocated(g%scal)) deallocate(g%scal)
 if (allocated(g%nres)) deallocate(g%nres)
 if (allocated(g%cres)) deallocate(g%cres)
end subroutine deallocategrid
