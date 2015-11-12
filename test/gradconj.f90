module grad_conj
  implicit none

  integer, parameter :: wp = kind(1.0D0)

contains

  subroutine grad_conj_opt(A,X,B)
    implicit none

    real(wp), dimension(:,:), intent(in)                  :: A
    real(wp), dimension(:), intent(in)                    :: B
    real(wp), dimension(size(A,1)), intent(out)           :: X
    real(wp), dimension(:,:), allocatable                 :: C
    real(wp), dimension(:), allocatable                   :: w, r, Aw
    real(wp)                                              :: alpha, beta, eps, err
    integer                                               :: i, j, n

    n = size(A,1)

    allocate(C(n,n), w(n), r(n), Aw(n))

    do i = 1, n
       C(i,i) = 1._wp !/ A(i,i)
    enddo

    eps = 0.0001_wp
    err = 1._wp

    X = 1._wp
    r = B - matmul(A,X)
    err = sqrt(dot_product(r,r))
    w = r
    alpha = dot_product(w,r) / dot_product(matmul(A,w),r)
    i = 0

    do while(eps<err)

       i = i + 1
       print*,i

       X = X + alpha*w
       r = B - matmul(A,X)
       Aw = matmul(A,w)
       err = sqrt(dot_product(r,r))
       beta = dot_product(Aw,r) / dot_product(Aw,w)
       w = r - beta*w
       alpha = dot_product(w,r) / dot_product(matmul(A,w),w)

       if(i==1000) then
          err = 0._wp
       endif

    enddo

  end subroutine grad_conj_opt

end module grad_conj


