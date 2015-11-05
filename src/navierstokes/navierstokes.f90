module navierstokes
     implicit none

     contains

!méthode de projection de Chorin
     subroutine projection_method(u,p,u_next,p_next,rho,nu,g,dt,dx,dy)
        implicit none

        ! vitesse au temps n
        real(kind=8) , dimension(:,:,2) , intent(in) :: u
        real(kind=8) , dimension(:,:) , intent(in) :: p
        
        !vitesse au temps n+1
        real(kind=8) , dimension(:,:,2), intent(out) :: u_next
        real(kind=8) , dimension(:,:) , intent(out) :: p_next

        !vitesse intermédiaire
        real(kind=8) , dimension(:,:) :: u_star

        !paramètres physiques
        real(kind=8) , intent(in) :: rho , nu , g

        !paramètres numériques
        real(kind=8) , intent(in) :: dt , dx , dy
        integer :: Nx , Ny

        real(kind=8) , dimension(:,:,2) , allocatable :: laplace_u , u_grad_u

        Nx = size(u,1)-1
        Ny = size(u,2)-1

        allocate(laplace_u(Nx-1,Ny-1,2),u_grad_u(Nx-1,Ny-1,2))

        laplace_u = ( u(3:Nx+1,2:Ny,:) - 2*u(2:Nx,2:Ny,:) + u(1:Nx-1,2:Ny,:) ) / (dx*dx) &
                   & + ( u(2:Nx,3:Ny+1,:) - 2*u(2:Nx,2:Ny,:) + u(2:Nx,1:Ny-1,:) ) / (dy*dy)
        
        u_grad_u = u(2:Nx,2:Ny,1)*( u(3:Nx+1,2:Ny,:) - u(1:Nx-1,2:Ny,:) ) / (2*dx) &
                  & + u(2:Nx,2:Ny,2)*( u(2:Nx,3:Ny+1,:) - u(2:Nx,1:Ny-1,:) ) / (2*dy)

        u_star = 0

        u_star(2:Nx,2:Ny,:) = u(2:Nx,2:Ny,:) + dt*( g + nu*laplace_u - u_grad_u )

        
        deallocate(laplace_u)

     end subroutine projection_method

!schéma MAC
     subroutine MAC_scheme()
        implicit none


     end subroutine MAC_scheme

!schéma Adams-Bashforth + Crank-Nicolson
     subroutine troisieme_schema()
        implicit none


     end subroutine troisieme_schema

module navierstokes
