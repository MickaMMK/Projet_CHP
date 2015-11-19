program test
  use grad_conj
  implicit none

  real(wp), dimension(1000,1000)    :: A
  real(wp), dimension(1000)         :: B, X
  integer                           :: i, j, n

  A = 0._wp

  n = 1000

  do i = 1, n
     A(i,i) = 4._wp
  enddo

  do i = 1, n-1
     A(i+1,i) = -1._wp
     A(i,i+1) = -1._wp
  enddo

  do i = 1, n-2
     A(i+2,i) = -1._wp
     A(i,i+2) = -1._wp
  enddo

  B = 1._wp

  call grad_conj_opt(A,X,B)

  X = matmul(A,X)
  X = X - B
  X = abs(X)
  print*,sum(X),'cest zéro ?'

end program test

