!
! Tioga is a library for overset grid assembly on parallel distributed systems
! Copyright (C) 2015 Jay Sitaraman
!
!=====================================================================
!
! read the grid files
!
!=====================================================================
subroutine readGrid_cell(g,myid)
!
use gridtype
implicit none
!
type(grid), intent(inout) :: g
integer, intent(in) :: myid
!
! local variables
!
character *128 :: fname,integer_string
character*1 :: hash
integer :: i,junk,j,m
real*8 :: rho,u,v,w,p,e,xx,yy,zz,rhoinf,pinf,uinf
real*8 :: fac,rsq
real*8 :: SIXTH=1./6
real*8 :: EIGHTH=1./8
!
write(integer_string,"(I7)") myid
fname='grid/cell'//trim(adjustl(integer_string))//'.plt'
!
open(unit=2,file=fname,form='formatted')
!
g%n4=0
g%n5=0
!
read(2,*) hash,g%n6,g%n8,g%nv,g%nCells,g%nwbc,g%nobc
!
read(2,*)
read(2,*)
read(2,*)
!
g%nvar=5
g%nghost=0
g%ndof=g%ncells+g%nghost
!
allocate(g%x(3*g%nv),g%bodytag(g%nv),g%iblank(g%nv))
allocate(g%scal(g%nvar))
allocate(g%xcentroid(3*g%ndof))
allocate(g%q(g%nvar*g%nv))
allocate(g%dq(3*g%nvar*g%nv))
!
g%iblank=1
!
! read coordinates and body tag 
! also create some conservative variables
!
rhoinf=1.
uinf=0.5
pinf=1/1.4
!
do i=1,g%nv
   read(2,*) g%x(3*i-2),g%x(3*i-1),g%x(3*i),g%bodytag(i)   
   xx=g%x(3*i-2)
   yy=g%x(3*i-1)
   zz=g%x(3*i)
   rsq=xx**2+yy**2+zz**2
   fac=0.1*exp(-rsq*0.5)

   rho=rhoinf*(1-fac)
   u=uinf+fac*xx/(rsq+0.01)
   v=fac*yy/(rsq+0.01)
   w=fac*zz/(rsq+0.01)
   p=pinf*(1-fac*fac)
   e=p/0.4+0.5*rho*(u**2+v**2+w**2)

   g%q(g%nvar*(i-1)+1)=rho
   g%q(g%nvar*(i-1)+2)=rho*u
   g%q(g%nvar*(i-1)+3)=rho*v
   g%q(g%nvar*(i-1)+4)=rho*w
   g%q(g%nvar*(i-1)+5)=e
enddo
!
g%scal=1
!
allocate(g%ndc6(6,g%n6),g%ndc8(8,g%n8))
!
g%ndc6=0
g%ndc8=0
g%nmax=8
!
! read prizm connectivity
!
m=0
do i=1,g%n6
 read(2,*) g%ndc6(1,i),g%ndc6(2,i),g%ndc6(3,i),junk,&
           g%ndc6(4,i),g%ndc6(5,i),g%ndc6(6,i),junk
 m=m+1
 g%xcentroid(3*m-2:3*m)=0.
 do j=1,6
    g%xcentroid(3*m-2)=g%xcentroid(3*m-2)+g%x(3*g%ndc6(j,i)-2)*SIXTH
    g%xcentroid(3*m-1)=g%xcentroid(3*m-1)+g%x(3*g%ndc6(j,i)-1)*SIXTH
    g%xcentroid(3*m)=g%xcentroid(3*m)+g%x(3*g%ndc6(j,i))*SIXTH
 enddo
enddo
!
! read hex connectivity
!
do i=1,g%n8
   read(2,*) g%ndc8(1,i),g%ndc8(2,i),g%ndc8(3,i),g%ndc8(4,i),&
             g%ndc8(5,i),g%ndc8(6,i),g%ndc8(7,i),g%ndc8(8,i)
    m=m+1
    g%xcentroid(3*m-2:3*m)=0.
    do j=1,8
       g%xcentroid(3*m-2)=g%xcentroid(3*m-2)+g%x(3*g%ndc8(j,i)-2)*EIGHTH
       g%xcentroid(3*m-1)=g%xcentroid(3*m-1)+g%x(3*g%ndc8(j,i)-1)*EIGHTH
       g%xcentroid(3*m)=g%xcentroid(3*m)+g%x(3*g%ndc8(j,i))*EIGHTH
    enddo
enddo
!
! read wall and overset boundary nodes
!
allocate(g%wbcnode(g%nwbc),g%obcnode(g%nobc))
!
do i=1,g%nwbc
   read(2,*) g%wbcnode(i)
enddo
!
do i=1,g%nobc
   read(2,*) g%obcnode(i)
enddo
!
return
!
end subroutine readGrid_cell
!=====================================================================
!
! read the grid files
!
!=====================================================================
subroutine readGrid_from_file(fname,g)
!
use gridtype
implicit none
!
character*(*) :: fname
type(grid), intent(inout) :: g
!
! local variables
!
character*128 :: a
character*1 :: hash
integer :: i,junk,j,m,k
real*8 :: rho,u,v,w,p,e,xx,yy,zz,rhoinf,pinf,uinf,qj
real*8 :: fac,rsq
real*8 :: SIXTH=1./6
real*8 :: EIGHTH=1./8
!
open(unit=101,file=fname,form='formatted')
!
do i=1,3
   read(101,"(A128)") a
enddo
do j=1,127
   if (a(j:j+1).eq.'N=' ) then
      k=j+2
      do while(a(k:k).ne.' ')
         k=k+1
      enddo
      read(a(j+2:k),*) g%nv
   endif
   if (a(j:j+1).eq.'E=') then
      k=j+2
     do while(a(k:k).ne.' ')
        k=k+1
     enddo
     read(a(j+2:k),*) g%ncells
   endif
enddo
g%n6=0
g%n8=g%nCells
g%nvar=1
g%nghost=0
g%ndof=g%ncells+g%nghost
!
allocate(g%x(3*g%nv),g%bodytag(g%nv),g%iblank(g%nv))
allocate(g%scal(g%nvar))
allocate(g%xcentroid(3*g%ndof))
allocate(g%q(g%nvar*g%nv))
allocate(g%dq(3*g%nvar*g%nv))
!
g%iblank=1
!
! read coordinates and body tag 
! also create some conservative variables
!
rhoinf=1.
uinf=0.5
pinf=1/1.4
!
do i=1,g%nv
   read(101,*) g%x(3*i-2),g%x(3*i-1),g%x(3*i),junk,g%bodytag(i) !,qj
enddo
!
g%scal=1
!
allocate(g%ndc6(6,g%n6),g%ndc8(8,g%n8))
!
g%ndc6=0
g%ndc8=0
g%nmax=8
!
! read prizm connectivity
!
m=0
do i=1,g%n6
 read(101,*) g%ndc6(1,i),g%ndc6(2,i),g%ndc6(3,i),junk,&
           g%ndc6(4,i),g%ndc6(5,i),g%ndc6(6,i),junk
 m=m+1
 g%xcentroid(3*m-2:3*m)=0.
 do j=1,6
    g%xcentroid(3*m-2)=g%xcentroid(3*m-2)+g%x(3*g%ndc6(j,i)-2)*SIXTH
    g%xcentroid(3*m-1)=g%xcentroid(3*m-1)+g%x(3*g%ndc6(j,i)-1)*SIXTH
    g%xcentroid(3*m)=g%xcentroid(3*m)+g%x(3*g%ndc6(j,i))*SIXTH
 enddo
enddo
!
! read hex connectivity
!
do i=1,g%n8
   read(101,*) g%ndc8(1,i),g%ndc8(2,i),g%ndc8(3,i),g%ndc8(4,i),&
             g%ndc8(5,i),g%ndc8(6,i),g%ndc8(7,i),g%ndc8(8,i)
    m=m+1
    g%xcentroid(3*m-2:3*m)=0.
    do j=1,8
       g%xcentroid(3*m-2)=g%xcentroid(3*m-2)+g%x(3*g%ndc8(j,i)-2)*EIGHTH
       g%xcentroid(3*m-1)=g%xcentroid(3*m-1)+g%x(3*g%ndc8(j,i)-1)*EIGHTH
       g%xcentroid(3*m)=g%xcentroid(3*m)+g%x(3*g%ndc8(j,i))*EIGHTH
    enddo
enddo
!
! read wall and overset boundary nodes
!
read(101,*) g%nwbc
allocate(g%wbcnode(g%nwbc))
!
do i=1,g%nwbc
   read(101,*) g%wbcnode(i)
enddo
read(101,*) g%nobc
allocate(g%obcnode(g%nobc))
!
do i=1,g%nobc
   read(101,*) g%obcnode(i)
enddo
!
close(101)
!
return
!
end subroutine readGrid_from_file



!=====================================================================
!
! move the grid
!
!=====================================================================
subroutine moveGrid(g)
use gridtype
implicit none
!
type(grid) :: g
!
integer :: i
real*8, save :: t=0
real*8, save :: period=4.
!
! begin
!
t=t+2*acos(-1.)/period
do i=1,g%nv
   if (g%bodytag(i)==1) then
      !
      ! translate the grid along the z-coordinate
      !
      g%x(3*i)=g%x(3*i)+0.04 !0.3*sin(t) 
   endif
enddo
!
do i=1,g%n6
  g%xcentroid(3*i)=g%xcentroid(3*i)+0.04
enddo
!
return
end subroutine moveGrid
