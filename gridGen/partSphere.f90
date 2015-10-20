!
! Tioga is a library for overset grid assembly on parallel distributed systems
! Copyright (C) 2015 Jay Sitaraman
!
subroutine getSpherePart(myid,numprocs,arange)
!
integer, intent(in) :: myid
integer, intent(in) :: numprocs
real*8, intent(out) :: arange(4)
!
integer :: i,nc,j,k
real*8 :: Ap,thetap,delta,Ai,a,thpolar
real*8, allocatable :: y(:)
integer, allocatable :: m(:)
real*8 :: pi
!
! begin
!
pi=acos(-1.)
Ap=4*pi/numprocs
!
thetap=asin(1-2.0/numprocs)
thpolar=pi*0.5-thetap
!
nc=nint((pi-2*thpolar)/sqrt(Ap))
delta=(pi-2*thpolar)/nc
!
allocate(y(nc),m(nc))
!
do i=1,nc
   Ai=2*pi*(sin(thetap-(i-1)*delta)-sin(thetap-i*delta))
   y(i)=Ai/Ap
enddo
!
a=0
do i=1,nc
   m(i)=nint(y(i)+a)
   a=a+(y(i)-m(i))
enddo
!
k=1
if (myid==k) then
  arange(1)=thetap
  arange(2)=pi*0.5
  arange(3)=0.
  arange(4)=2*pi
endif
!
do i=1,nc
   do j=1,m(i)
     k=k+1
     if (myid==k) then
        arange(1)=thetap-i*delta
	arange(2)=arange(1)+delta
        arange(3)=(j-1)*2*pi/m(i)
        arange(4)=arange(3)+2*pi/m(i)
	return
     endif
   enddo
enddo
!
k=numprocs
if (myid==k) then
  arange(1)=-pi*0.5
  arange(2)=-thetap
  arange(3)=0.
  arange(4)=2*pi
endif
!
!
return
end subroutine getSpherePart
