module grad_conj
  implicit none

  integer, parameter :: wp = kind(1.0D0)

contains

  subroutine grad_conj_opt(X,B,dx)
    implicit none

!    real(wp), dimension(:,:), intent(in)                  :: A
    real(wp), dimension(:), intent(in)                    :: B
    real(wp), dimension(size(B,1)), intent(out)           :: X
    real(wp), dimension(:,:), allocatable                 :: C, A_modif
    real(wp), dimension(:), allocatable                   :: w, r, Aw, B_modif
    real(wp), intent(in)                                  :: dx
    real(wp)                                              :: alpha, beta, eps, err
    integer                                               :: i, j, n

!    n = size(A,1)

   ! allocate(C(n,n), w(n), r(n), Aw(n))
    allocate(w(n), r(n), Aw(n))

!!$    C = 0._wp
!!$
!!$
!!$    do i = 1, n
!!$       C(i,i) = 1._wp / A(i,i) ! matmul(A,B) = mat_vect(B)
!!$    enddo
!!$
!!$    A_modif = matmul(C,A)
!!$    B_modif = matmul(C,B)

    eps = 0.0000001_wp
    err = 1._wp

    X = 1._wp
    r = B - mat_vect(X)
    err = sqrt(dot_product(r,r))
    w = r
    alpha = dot_product(w,r) / dot_product(mat_vect(w),r)
    i = 0

    do while((eps<err).and.(i<10000))

       i = i + 1
       X = X + alpha*w
       r = B - mat_vect(X)
       Aw = mat_vect(w)
       err = sqrt(dot_product(r,r))
       beta = dot_product(Aw-w,r) / dot_product(Aw,w)
       w = r - beta*w
       alpha = dot_product(w,r) / dot_product(mat_vect(w),w)

       print*,err

    enddo

    ! X = X * dx * dx

    print*,i

  end subroutine grad_conj_opt








  function mat_vect(X) result(Y)
    implicit none

    real(kind=8) , dimension(:) , intent(in)   :: X
    real(kind=8) , dimension(size(X))          :: Y

    integer :: N , i , j

    N = int(sqrt(float(size(X))))

    Y(1) = -2*X(1) + X(2) + X(N+1)
    Y(1) = -Y(1)/2

    do j = 2 , N-1

       Y(j) = -3*X(j) + X(j-1) + X(j+1) + X(j+N)
       Y(j) =-Y(j)/3

    end do

    Y(N) = -2*X(N) + X(N-1) + X(2*N)
    Y(N) = -Y(N)/2

    do i = 2 , N-1

       Y((i-1)*N+1) = -3*X((i-1)*N+1) + X((i-1)*N+2) + X((i-1)*N+N+1) + X((i-1)*N-N+1) 
       Y((i-1)*N+1) =-Y((i-1)*N+1)/3

       do j = 2 , N-1

          Y((i-1)*N+j) = -4*X((i-1)*N+j) + X((i-1)*N+j+1) + X((i-1)*N+j-1) &
               &+ X((i-1)*N+j-N) + X((i-1)*N+j+N)
          Y((i-1)*N+j) =-Y((i-1)*N+j)/4

       end do

       Y(i*N) = -3*X(i*N) + X(i*N-1) + X((i-1)*N) + X((i+1)*N) 
       Y(i*N) =-Y(i*N)/3

    end do

    Y(N*(N-1)+1) = -2*X(N*(N-1)+1) + X(N*(N-1) + 2) + X(N*(N-1)-N+1)
    Y(N*(N-1)+1)= -Y(N*(N-1)+1)/2

    do j = 2 , N-1

       Y(N*(N-1)+j) = -3*X(N*(N-1)+j) + X(N*(N-1)+j+1) + X(N*(N-1)+j-1)+ X(N*(N-1)+j-N)
       Y(N*(N-1)+j) =-Y(N*(N-1)+j)/3

    end do

    ! on remplace la dernière ligne de la matrice par (0--------01) pour fixer la valeur de la pression au centre (N,N)
    !Y(N*N) = -2*X(N*N) + X(N*N-1) + X(N*N-N)              
    Y(N*N) = X(N*N)

    !  deallocate(Y)

  end function mat_vect











!!$subroutine grad_conj_test(A,X,B)
!!$    implicit none
!!$
!!$   real(wp), dimension(:,:), intent(in)                  :: A
!!$    real(wp), dimension(:), intent(in)                    :: B
!!$    real(wp), dimension(size(A,1)), intent(out)           :: X
!!$    real(wp), dimension(:,:), allocatable                 :: Id
!!$    real(wp), dimension(:), allocatable                   :: g1, g2, d
!!$    real(wp)                                              :: alpha, beta, eps, err, w1
!!$    integer                                               :: i, j, n
!!$
!!$    n = size(A,1)
!!$
!!$    allocate(Id(n,n), g1(n), g2(n), d(n))
!!$
!!$    Id = 0._wp
!!$    do i = 1, n
!!$       Id(i,i) = 1._wp
!!$    enddo
!!$
!!$    w1 = 0.95_wp
!!$
!!$
!!$    X = -B
!!$    g1 = matmul(A,X) - B
!!$    d = g1
!!$
!!$    err = 1.0_wp
!!$
!!$    
!!$    do while(err>0.0001_wp)
!!$       alpha = 0.001_wp
!!$    !   do while(dot_product(d,(matmul((1._wp-w1)*A+Id),X)+alpha*d+(w1-1._wp)*B)>0._wp)
!!$    !      alpha = alpha / 2._wp
!!$    !   enddo
!!$
!!$       X = X + alpha*d
!!$       g2 = matmul(A,X) - B
!!$       beta = dot_product(g2,g2-g1) / dot_product(g1,g1)
!!$       d = -g2 + beta*d
!!$       g1 = g2
!!$
!!$       err = sqrt(dot_product(g2,g2))
!!$
!!$       print*,err
!!$
!!$    enddo
!!$
!!$
!!$
!!$  end subroutine grad_conj_test




end module grad_conj


