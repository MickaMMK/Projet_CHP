module navierstokes
     implicit none

     contains

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
     subroutine projection_method(u,p,u_next,p_next,rho,nu,g,dt,dx,dy)
        implicit none

        ! vitesse au temps n
        real(kind=8) , dimension(:,:,:) , intent(in) :: u
        real(kind=8) , dimension(:,:) , intent(in) :: p
        
        !vitesse au temps n+1
        real(kind=8) , dimension(:,:,:), intent(out) :: u_next
        real(kind=8) , dimension(:,:) , intent(out) :: p_next

        !vitesse intermédiaire
        real(kind=8) , dimension(:,:,:) , allocatable :: u_star

        !paramètres physiques
        real(kind=8) , intent(in) :: rho , nu , g

        !paramètres numériques
        real(kind=8) , intent(in) :: dt , dx , dy
        integer :: Nx , Ny

        real(kind=8) , dimension(:,:,:) , allocatable :: laplace_u , u_grad_u

        Nx = size(u,1)-1
        Ny = size(u,2)-1

        allocate(laplace_u(Nx-1,Ny-1,2),u_grad_u(Nx-1,Ny-1,2))
        allocate(u_star(Nx,Ny,2))

        laplace_u = ( u(3:Nx+1,2:Ny,:) - 2*u(2:Nx,2:Ny,:) + u(1:Nx-1,2:Ny,:) ) / (dx*dx) &
                   & + ( u(2:Nx,3:Ny+1,:) - 2*u(2:Nx,2:Ny,:) + u(2:Nx,1:Ny-1,:) ) / (dy*dy)
        
      !  u_grad_u = u(2:Nx,2:Ny,1)*( u(3:Nx+1,2:Ny,:) - u(1:Nx-1,2:Ny,:) ) / (2*dx) &
      !            & + u(2:Nx,2:Ny,2)*( u(2:Nx,3:Ny+1,:) - u(2:Nx,1:Ny-1,:) ) / (2*dy)

        u_star = 0

        u_star(2:Nx,2:Ny,:) = u(2:Nx,2:Ny,:) + dt*( g + nu*laplace_u - u_grad_u )

        !résolution problème de Poisson 

                !remplissage matrice


                !remplissage second membre


                !gradient conjugué



        !calcul vitesse au temps n+1
        
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
