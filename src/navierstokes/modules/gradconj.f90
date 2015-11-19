module grad_conj
  implicit none

  integer, parameter :: wp = kind(1.0D0)

contains

  subroutine grad_conj_opt(A,X,B)
    implicit none

   !real(wp), dimension(:,:), intent(in)                  :: A
    real(wp), dimension(:), intent(in)                    :: B
    real(wp), dimension(size(A,1)), intent(out)           :: X
    real(wp), dimension(:,:), allocatable                 :: C, A_modif
    real(wp), dimension(:), allocatable                   :: w, r, Aw, B_modif
    real(wp)                                              :: alpha, beta, eps, err
    integer                                               :: i, j, n

    n = size(A,1)

    allocate(C(n,n), w(n), r(n), Aw(n))

!!$    do i = 1, n
!!$       C(i,i) = 100._wp / A(i,i) ! matmul(A,B) = mat_vect(B)
!!$    enddo
!!$
!!$    A_modif = matmul(C,A)
!!$    B_modif = matmul(C,B)

    eps = 0.0001_wp
    err = 1._wp

    X = 1._wp
    r = B - matmul(A,X)
    err = sqrt(dot_product(r,r))
    w = r
    alpha = dot_product(w,r) / dot_product(mat_vect(w),r)
    i = 0

    do while((eps<err).and.(i<1000))

       i = i + 1
       X = X + alpha*w
       r = B - mat_vect(X)
       Aw = mat_vect(w)
       err = sqrt(dot_product(r,r))
       beta = dot_product(Aw-w,r) / dot_product(Aw,w)
       w = r - beta*w
       alpha = dot_product(w,r) / dot_product(mat_vect(w),w)

    enddo

    print*,i

  end subroutine grad_conj_opt



subroutine grad_conj_test(A,X,B)
    implicit none

   !real(wp), dimension(:,:), intent(in)                  :: A
    real(wp), dimension(:), intent(in)                    :: B
    real(wp), dimension(size(A,1)), intent(out)           :: X
    real(wp), dimension(:,:), allocatable                 :: C, A_modif
    real(wp), dimension(:), allocatable                   :: w, r, Aw, B_modif
    real(wp)                                              :: alpha, beta, eps, err
    integer                                               :: i, j, n

    n = size(A,1)

   ! X = 



  end subroutine grad_conj_test




end module grad_conj


