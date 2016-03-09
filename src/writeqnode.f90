subroutine writeqnode(myid,qnode,qnodesize)
implicit none
integer :: i
integer, intent(in) :: myid,qnodesize
real*8, intent(in) :: qnode(qnodesize)

write(1000+myid,1000) (qnode(i),i=1,qnodesize)
1000 format(200(1x,E14.8))
return
end subroutine writeqnode
