module grad_conj
  implicit none

  integer, parameter :: wp = kind(1.0D0)

contains

  subroutine grad_conj_opt(A,X,B)
    implicit none

    real(wp), dimension(:,:), intent(in)                  :: A
    real(wp), dimension(:), intent(in)                    :: B
    real(wp), dimension(size(A,1)), intent(out)           :: X
    real(wp), dimension(:,:), allocatable                 :: C, A_modif
    real(wp), dimension(:), allocatable                   :: w, r, Aw, B_modif
    real(wp)                                              :: alpha, beta, eps, err
    integer                                               :: i, j, n

    n = size(A,1)

    allocate(C(n,n), w(n), r(n), Aw(n))

    do i = 1, n
       C(i,i) = 100._wp / A(i,i) ! matmul(A,mathieu) = produit_mat_vect(X)
    enddo

    A_modif = matmul(C,A)
    B_modif = matmul(C,B)

    eps = 0.0001_wp
    err = 1._wp

    X = 1._wp
    r = B_modif - matmul(A_modif,X)
    err = sqrt(dot_product(r,r))
    w = r
    alpha = dot_product(w,r) / dot_product(matmul(A_modif,w),r)
    i = 0

    do while((eps<err).and.(i<1000))

       i = i + 1
       X = X + alpha*w
       r = B_modif - matmul(A_modif,X)
       Aw = matmul(A_modif,w)
       err = sqrt(dot_product(r,r))
       beta = dot_product(Aw-w,r) / dot_product(Aw,w)
       w = r - beta*w
       alpha = dot_product(w,r) / dot_product(matmul(A_modif,w),w)

    enddo

    print*,i

  end subroutine grad_conj_opt

end module grad_conj


