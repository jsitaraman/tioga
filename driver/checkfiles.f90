!
! Tioga is a library for overset grid assembly on parallel distributed systems
! Copyright (C) 2015 Jay Sitaraman
!
program checkfiles
 integer, parameter :: n = 10340
 integer, parameter :: p = 4
 integer a(p,n)
 integer b(p,n)
 real*8  da(p,n)
 real*8  db(p,n)

 read(200,*) a
 read(201,*) da
 open(unit=1,file='adtInts.dat',form='formatted')
 read(1,*) b
 close(1)
 open(unit=1,file='adtReals.dat',form='formatted')
 read(1,*) db
 close(1)

 isum=0
 dsum=0
 do i=1,n
  do j=1,p
   isum=isum+abs(a(j,i)-b(j,i)-1)
   dsum=dsum+abs(da(j,i)-db(j,i))
  enddo
 enddo

 write(6,*) isum
 write(6,*) dsum
end 
