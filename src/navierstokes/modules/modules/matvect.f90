module matvectmod

  implicit none

contains

  function mat_vect(X,dx) result(Y)
    implicit none

    real(kind=8) , dimension(:) , intent(in)   :: X
    real(kind=8) , intent(in)                  :: dx
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
    Y(N*N) = X(N*N)*dx

    !  deallocate(Y)

  end function mat_vect

end module matvectmod
