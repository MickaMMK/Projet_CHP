module navierstokes
  use grad_conj
  implicit none

contains

  !produit matrice vecteur 


  function mat_vect(X,dx) result(Y)
    implicit none

    real(kind=8) , dimension(:) , intent(in) :: X
    real(kind=8) :: dx
    real(kind=8) , dimension(size(X))  :: Y

    integer :: N , i , j

    N = int(sqrt(float(size(X))))

    Y(1) = -2*X(1) + X(2) + X(N+1)

    do j = 2 , N-1

       Y(j) = -3*X(j) + X(j-1) + X(j+1) + X(j+N)

    end do

    Y(N) = -2*X(N) + X(N-1) + X(2*N)

    do i = 2 , N-1

       Y((i-1)*N+1) = -3*X((i-1)*N+1) + X((i-1)*N+2) + X((i-1)*N+N+1) + X((i-1)*N-N+1) 

       do j = 2 , N-1

          Y((i-1)*N+j) = -4*X((i-1)*N+j) + X((i-1)*N+j+1) + X((i-1)*N+j-1) &
               &+ X((i-1)*N+j-N) + X((i-1)*N+j+N)

       end do

       Y(i*N) = -3*X(i*N) + X(i*N-1) + X(i*N+1) + X((i+1)*N) 

    end do

    Y(N*(N-1)+1) = -2*X(N*(N-1)+1) + X(N*(N-1) + 2) + X(N*(N-1)-N+1)

    do j = 2 , N-1

       Y(N*(N-1)+j) = -3*X(N*(N-1)+j) + X(N*(N-1)+j+1) + X(N*(N-1)+j-1)+ X(N*(N-1)+j-N)

    end do

    ! on remplace la dernière ligne de la matrice par (0--------01) pour fixer la valeur de la pression au centre (N,N)
    !Y(N*N) = -2*X(N*N) + X(N*N-1) + X(N*N-N)              
    Y(N*N) = X(N*N)

    Y = (1/(dx*dx))*Y

    !  deallocate(Y)

  end function mat_vect

  !remplissage de la matrice du problème de Poisson
  subroutine remplissage_poisson(A,dx,N)
    implicit none

    real(kind=8) , dimension(:,:) , intent(out) :: A
    real(kind=8) , intent(in) :: dx
    integer , intent(in) :: N !matrice de pression de taille NxN

    integer :: i , j 

    A = 0

    !premier bloc de N lignes

    A(1,1) = -2
    A(1,2) = 1
    A(1,N+1) = 1

    do i = 2 , N-1 

       A(i,i) = -3
       A(i,i+1) = 1
       A(i,i-1) = 1
       A(i,i+N) = 1

    end do

    A(N,N) = -2
    A(N,N-1) = 1
    A(N,2*N) = 1

    !dernier bloc de N lignes

    A(N*(N-1)+1,N*(N-1)+1) = -2
    A(N*(N-1)+1,N*(N-1)+2) = 1
    A(N*(N-1)+1,N*(N-1)+1-N) = 1

    do i = 2 , N-1 

       A(N*(N-1)+i,N*(N-1)+i) = -3
       A(N*(N-1)+i,N*(N-1)+i+1) = 1
       A(N*(N-1)+i,N*(N-1)+i-1) = 1
       A(N*(N-1)+i,N*(N-1)+i-N) = 1

    end do

    A(N*N,N*N) = -2
    A(N*N,N*N-1) = 1
    A(N*N,N*N-N) = 1

    !blocs génériques

    do i = 2 , N-1

       A((i-1)*N+1,(i-1)*N+1) = -3
       A((i-1)*N+1,(i-1)*N+2) = 1
       A((i-1)*N+1,(i-1)*N+1+N) = 1
       A((i-1)*N+1,(i-1)*N+1-N) = 1

       do j = 2 , N-1

          A((i-1)*N+j,(i-1)*N+j) = -4
          A((i-1)*N+j,(i-1)*N+j+1) = 1
          A((i-1)*N+j,(i-1)*N+j-1) = 1
          A((i-1)*N+j,(i-1)*N+j+N) = 1
          A((i-1)*N+j,(i-1)*N+j-N) = 1

       end do

       A(i*N,i*N) = -3
       A(i*N,i*N-1) = 1
       A(i*N,i*N+N) = 1
       A(i*N,i*N-N) = 1

    end do

    A = (1./(dx*dx))*A

  end subroutine remplissage_poisson

  !méthode de projection de Chorin
  subroutine projection_method(u,p,rho,rho_centre,nu,g,dt,dx)
    implicit none

    ! vitesse au temps n
    real(kind=8), dimension(:,:,:), intent(inout)       :: u
    real(kind=8), dimension(:,:), intent(inout)         :: p

    !vitesse au temps n+1
    real(kind=8), dimension(size(u,1),size(u,2),size(u,3))                   :: u_next
    real(kind=8), dimension(size(p,1),size(p,2))                             :: p_next

    !vitesse intermédiaire
    real(kind=8), dimension(:,:,:), allocatable      :: u_star

    !paramètres physiques
    real(kind=8), intent(in)                         :: nu, g
    real(kind=8), dimension(:,:), intent(in)         :: rho, rho_centre

    !paramètres numériques
    real(kind=8), intent(in)                         :: dt, dx 
    integer :: N, i, j

    !matrice de Poisson et second membre
    real(kind=8), dimension(:), allocatable          :: B
    real(kind=8), dimension(:,:), allocatable        :: A

    real(kind=8), dimension(:,:,:), allocatable      :: laplace_u , u_grad_u , grad_p

    N = size(u,1)-1

    allocate(laplace_u(2:N,2:N,2),u_grad_u(2:N,2:N,2),grad_p(2:N,2:N,2))
    allocate(u_star(N+1,N+1,2))
    allocate(B(N))
    allocate(A(N,N))

    laplace_u = ( u(3:N+1,2:N,:) - 2*u(2:N,2:N,:) + u(1:N-1,2:N,:) ) / (dx*dx) &
         & + ( u(2:N,3:N+1,:) - 2*u(2:N,2:N,:) + u(2:N,1:N-1,:) ) / (dx*dx)

    do i = 2 , N
       do j = 2 , N

          u_grad_u(i,j,:) = u(i,j,1)*( u(i+1,j,:) - u(i-1,j,:) ) / (2*dx) &
               & + u(i,j,2)*( u(i,j+1,:) - u(i,j-1,:) ) / (2*dx)         

       end do
    end do

    u_star = 0

    u_star(2:N,2:N,:) = u(2:N,2:N,:) + dt*( g + nu(2:N,2:N)*laplace_u - u_grad_u )

    !résolution problème de Poisson 

    !remplissage matrice
    call remplissage_poisson(A,dx,N)

    !remplissage second membre
    do i = 1 , N
       do j = 1 , N

          B(i+N*(j-1)) = (rho_centre(i,j)/(2*dt*dx))*(u_star(i+1,j+1,1)+u_star(i+1,j,1)-u_star(i,j+1,1)-u_star(i,j,1) &
               & +u_star(i,j+1,2)+u_star(i+1,j+1,2)-u_star(i+1,j,2)-u_star(i,j,2))

       end do
    end do
    !gradient conjugué

    call grad_conj_opt(P_next,B,dx)

    !calcul vitesse au temps n+1

    do i = 2 , N

       do j = 2 , N

          grad_p(i,j,1) =  p_next(i,j)+p_next(i,j-1)-p_next(i-1,j)-p_next(i-1,j-1)
          grad_p(i,j,2) =  p_next(i,j)+p_next(i-1,j)-p_next(i-1,j-1)-p_next(i,j-1)

       end do

    end do

    u_next = 0

    u_next(2:N,2:N,:) = u_star(2:N,2:N,:) - (dt/(2*dx*rho(2:N,2:N)))*grad_p(2:N,2:N,:)  

    u = u_next
    p = p_next
    
    deallocate(laplace_u,u_grad_u)

  end subroutine projection_method

  !schéma MAC
  subroutine MAC_scheme()
    implicit none


  end subroutine MAC_scheme

  !schéma Adams-Bashforth + Crank-Nicolson
  subroutine troisieme_schema()
    implicit none


  end subroutine troisieme_schema

end module navierstokes
