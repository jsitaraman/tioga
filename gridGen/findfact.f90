!
! Tioga is a library for overset grid assembly on parallel distributed systems
! Copyright (C) 2015 Jay Sitaraman
!
subroutine factorize(n,ff)
  implicit none
  integer, intent(in) :: n
  integer, intent(out) :: ff(3)
  integer, allocatable, dimension(:) :: factors
  integer :: i, f,k
  real*8 :: nc
  !
  if (n==1) then
   ff=1
   return
  endif
  allocate(factors(n/2))
  !
  call primefactors(n,factors,f)
  !
  if (f .ge. 3) then
     nc=(1.0*n)**(1./3.)
     k=1 
     ff(k)=1
     i=1
     do while(i.le.f .and. k.le.2)
        ff(k)=ff(k)*factors(i)
        if (ff(k) > 1.8*nc) then
           ff(k)=ff(k)/factors(i)
           k=k+1
           ff(k)=1
           i=i-1
        endif
        i=i+1
     enddo
     ff(3)=1
     do while(i.le.f) 
        ff(3)=ff(3)*factors(i)
        i=i+1
     enddo
  elseif (f.eq.2) then
     ff(1)=factors(1)
     ff(2)=factors(2)
     ff(3)=1
  else
     ff(1)=factors(1)
     ff(2)=1
     ff(3)=1
  endif
  deallocate(factors)
end subroutine factorize
!
!subroutine to find the prime factors of a number
!
subroutine primefactors(num, factors, f)
!
  implicit none
  integer, intent(in) :: num  !input number
  integer,intent(out), dimension((num/2))::factors !array to store factors
  integer, intent(inout) :: f
  integer :: i, n
  i = 2  !eligible factor
  f = 1  !number of factors
  n = num !store input number into a temporary variable
  do
     if (mod(n,i) == 0) then !if i divides 2, it is a factor
        factors(f) = i
        f = f+1
        n = n/i
     else
        i = i+1     !not a factor. move to next number
     end if
     if (n == 1) then		
        !since f is incremented after a factor is found
        f = f-1		!its value will be one more than the number of factors
        !hence the value of f is decremented
        exit
     end if
  end do
end subroutine primefactors
