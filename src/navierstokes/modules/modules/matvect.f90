module matvectmod

  implicit none

contains

  ! fonction qui à X associe AX avec A la matrice du problème de Poisson en pression en monophasique

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


  ! fonction qui à X associe AX avec A matrice du problème de Poisson en pression en diphasique

  function mat_vect_diphasique(X,dx,rho) result(Y)
    implicit none

    real*8, dimension(:), intent(in)              :: X
    real*8, intent(in)                            :: dx
    real*8, dimension(:), intent(in)              :: rho

    real*8, dimension(:), allocatable             :: Y

    integer                                       :: i, j, N

    N = int(sqrt(float(size(X))))

    allocate(Y(size(X)))

    !premier bloc de N lignes

    Y(1) = -X(1)*(1./(rho(1)+rho(2)) + 1./(rho(1)+rho(N+1))) + X(2)/(rho(2)+rho(1)) + X(N+1)/(rho(N+1)+rho(1))
    ! Y(1) = -Y(1)/2

    do j = 2 , N-1

       Y(j) = -X(j)*( 1./(rho(j)+rho(j-1)) + 1./(rho(j)+rho(j+1)) &
            & + 1./(rho(j)+rho(j+N)) ) &
            & + X(j-1)/(rho(j)+rho(j-1)) & 
            & + X(j+1)/(rho(j)+rho(j+1)) &
            & + X(j+N)/(rho(j)+rho(j+N))
       !   Y(j) =-Y(j)/3

    end do

    Y(N) = -X(N)*(1/(rho(N)+rho(N-1)) + 1/(rho(N)+rho(2*N))) + X(N-1)/(rho(N)+rho(N-1)) + X(2*N)/(rho(N)+rho(2*N)) 
    !Y(N) = -Y(N)/2


    ! N-2 blocs de N lignes

    do i = 2 , N-1

       Y((i-1)*N+1) = - X((i-1)*N+1) * (1./(rho((i-1)*N+1)+rho((i-1)*N+2)) &
            & + 1./(rho((i-1)*N+1)+rho((i-1)*N+N+1)) + 1./(rho((i-1)*N+1)+rho((i-1)*N-N+1))) &
            & + X((i-1)*N+2)/(rho((i-1)*N+1)+rho((i-1)*N+2)) &
            & + X((i-1)*N+N+1)/(rho((i-1)*N+1)+rho((i-1)*N+N+1)) &
            & + X((i-1)*N-N+1)/(rho((i-1)*N+1)+rho((i-1)*N-N+1))
       !  Y((i-1)*N+1) =-Y((i-1)*N+1)/3


       do j = 2 , N-1

          Y((i-1)*N+j) = -X((i-1)*N+j) * ( 1./(rho((i-1)*N+j  ) + rho((i-1)*N+j+1  )   )  &
               & + 1./( rho( (i-1)*N+j ) + rho((i-1)*N+j-1  )    ) &
               & + 1./(   rho( (i-1)*N+j ) + rho((i-1)*N+j-N  )      )  + 1./(  rho( (i-1)*N+j ) + rho((i-1)*N+j+N  )      )    ) &
               & + X((i-1)*N+j+1)/(rho((i-1)*N+j  ) + rho((i-1)*N+j+1  )) + X((i-1)*N+j-1)/(rho( (i-1)*N+j ) + rho((i-1)*N+j-1  )) &
               & + X((i-1)*N+j-N)/(rho( (i-1)*N+j ) + rho((i-1)*N+j-N  )) + X((i-1)*N+j+N)/(rho( (i-1)*N+j ) + rho((i-1)*N+j+N  ))
          !    Y((i-1)*N+j) =-Y((i-1)*N+j)/4

       end do

       Y(i*N) = -X(i*N) * (1./( rho( i*N )  + rho( i*N - 1  ) ) &
            & + 1./(  rho( i*N  )  + rho((i-1)*N  )  ) + 1./(  rho( i*N  ) + rho( (i+1)*N  )  )  ) &
            & + X(i*N-1)/( rho( i*N )  + rho( i*N - 1  )  ) &
            & + X((i-1)*N)/( rho( i*N  )  + rho((i-1)*N  )  ) &
            & + X((i+1)*N)/( rho( i*N  ) + rho( (i+1)*N  )  )   
       !Y(i*N) =-Y(i*N)/3

    end do


    ! dernier bloc de N lignes

    Y(N*(N-1)+1) = -X(N*(N-1)+1) * ( 1./( rho(N*(N-1)+1) + rho(N*(N-1)+2) ) &
         & + 1./(rho(N*(N-1)+1 ) + rho(N*(N-1)-N+1)) ) & 
         & + X(N*(N-1) + 2)/(rho(N*(N-1)+1) + rho(N*(N-1)+2 )) &
         & + X(N*(N-1)-N+1)/(rho(N*(N-1)+1) + rho(N*(N-1)-N+1)) 
    !Y(N*(N-1)+1)= -Y(N*(N-1)+1)/2

    do j = 2 , N-1

       Y(N*(N-1)+j) = -X(N*(N-1)+j) * ( 1./(rho( N*(N-1)+j) + rho(N*(N-1)+j+1)) &
            & + 1./(rho( N*(N-1)+j) + rho(N*(N-1)+j-1)) + 1./(  rho(N*(N-1)+j) + rho( N*(N-1)+j-N))) &
            & + X(N*(N-1)+j+1)/(rho( N*(N-1)+j) + rho(N*(N-1)+j+1)) &
            & + X(N*(N-1)+j-1)/(rho(N*(N-1)+j) + rho(N*(N-1)+j-1)) &
            & + X(N*(N-1)+j-N)/(rho(N*(N-1)+j) + rho( N*(N-1)+j-N))
       !  Y(N*(N-1)+j) =-Y(N*(N-1)+j)/3

    end do

    ! on remplace la dernière ligne de la matrice par (0--------01) pour fixer la valeur de la pression au centre (N,N)
    !Y(N*N) = -2*X(N*N) + X(N*N-1) + X(N*N-N)              
    !  Y(N*N) = X(N*N)*dx

    Y = 2 * Y / (dx*dx)    
    Y(N*N) = X(N*N)*dx

  end function mat_vect_diphasique

end module matvectmod
