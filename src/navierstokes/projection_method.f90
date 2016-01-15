!------------------------------------------------------------------------------
! ENSEIRB-MATMECA - Colonne Splash 2015(c) - All rights reserved
!------------------------------------------------------------------------------

module projection_methodmod
  use grad_conj
  use vectmod
  implicit none

  real(8), parameter :: Patm = 1.013d5

contains

   !---------------------------------------------------------------------------
   !> @brief
   !> Uses the projection method to compute \f$ u_{next} \f$ and \f$ p_{next} \f$.
   !
   !> @author 
   !> Corentin PRIGENT, student engineer at ENSEIRB-MATMECA, Bordeaux, FR.
   !
   !> @param[in]       
   !> @param[out]       
   !> @return 
   !--------------------------------------------------------------------------- 

  !méthode de projection de Chorin monophasique
  subroutine projection_method(u,p,rho,rho_centre,nu,g,dt,dx,level,A,ipvt)
    implicit none

    ! vitesse au temps n
    real(kind=8), dimension(:,:,:), intent(inout)       :: u
    real(kind=8), dimension(:,:), intent(inout)         :: p
    real(kind=8), dimension(:,:), intent(in)         :: level

    !vitesse au temps n+1
    real(kind=8), dimension(size(u,1),size(u,2),size(u,3))  :: u_next
    real(kind=8), dimension(size(p,1),size(p,2))            :: p_next
    real(kind=8), dimension(size(p,1)*size(p,2))            :: p_next_vect

    !vitesse intermédiaire
    real(kind=8), dimension(:,:,:), allocatable      :: u_star

    !paramètres physiques
    real(kind=8), dimension(2), intent(in)           :: g
    real(kind=8), dimension(:,:), intent(in)         :: rho, rho_centre, nu

    !paramètres numériques
    real(kind=8), intent(in)                         :: dt, dx 
    integer :: N, i, j

    !matrice de Poisson et second membre
    real(kind=8), dimension(:,:), intent(in)         :: A
    integer, dimension(:), intent(in)                :: ipvt
    real(kind=8), dimension(:), allocatable          :: B

    real(kind=8), dimension(:,:,:), allocatable      :: laplace_u , u_grad_u , grad_p

    integer                                          :: info, im, jm, km, k

    character(len=3)                                 :: nom

    real(8) :: data

    N = size(u,1)-1

    allocate(laplace_u(2:N,2:N,2),u_grad_u(2:N,2:N,2),grad_p(2:N,2:N,2))
    allocate(u_star(N+1,N+1,2))
    allocate(B(N*N))

    laplace_u = ( u(3:N+1,2:N,:) - 2*u(2:N,2:N,:) + u(1:N-1,2:N,:) ) / (dx*dx) &
         & + ( u(2:N,3:N+1,:) - 2*u(2:N,2:N,:) + u(2:N,1:N-1,:) ) / (dx*dx)

    do i = 2 , N
       do j = 2 , N

          u_grad_u(i,j,:) = u(i,j,1)*( u(i+1,j,:) - u(i-1,j,:) ) / (2*dx) &
               & + u(i,j,2)*( u(i,j+1,:) - u(i,j-1,:) ) / (2*dx)         

       end do
    end do

    u_star = 0
!!$    u_star(:,N+1,1) = 1.

    do i = 1, size(u_star,3)
       u_star(2:N,2:N,i) = u(2:N,2:N,i) + dt*( g(i) + nu(2:N,2:N)*laplace_u(:,:,i) - u_grad_u(:,:,i) )
    end do

    !résolution problème de Poisson 

!!$    call remplissage_poisson(A,dx,N)

    !remplissage second membre
    do i = 1 , N
       do j = 1 , N

          B(i+N*(j-1)) = (rho_centre(i,j)/(2*dt*dx))*(u_star(i+1,j+1,1)+u_star(i+1,j,1)-u_star(i,j+1,1)-u_star(i,j,1) &
               & +u_star(i,j+1,2)+u_star(i+1,j+1,2)-u_star(i+1,j,2)-u_star(i,j,2))

       end do
    end do


    !gradient conjugué

    B = condi(B*dx*dx)

    !!!!!!!!!!!!!!!!!!!!!!! CHOIX SOLVEUR !!!!!!!!!!!!!!!!!!!!!!!

    P_next_vect = vect1(P_next)
    call grad_conj_opt(P_next_vect,B,dx)
    P_next = unvect1(remoy(P_next_vect,Patm))
!!$
!!$    !!!call DGETRF(N*N, N*N, A, N*N, ipvt, info)      ! On n'utilise plus DGETRF (fait dans le main)
!!$    call DGETRS('N', N*N, 1, A, N*N, ipvt, B, N*N, info)
!!$    P_next = unvect1(remoy(B,Patm))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i = 2 , N

       do j = 2 , N

          grad_p(i,j,1) =  p_next(i,j)+p_next(i,j-1)-p_next(i-1,j)-p_next(i-1,j-1)
          grad_p(i,j,2) =  p_next(i,j)+p_next(i-1,j)-p_next(i-1,j-1)-p_next(i,j-1)

       end do

    end do

    u_next = 0
   ! u_next(:,N+1,1) = 0.

    do i = 1, size(u_next,3)
       u_next(2:N,2:N,i) = u_star(2:N,2:N,i) - (dt/(2*dx*rho(2:N,2:N)))*grad_p(2:N,2:N,i)
    end do

    im = 1
    jm = 1
    km = 1
    do i = 1, N
       do j = 1, N
          do k = 1, 2
             if(abs(u_next(i,j,k)) > abs(u_next(im,jm,km))) then
                im = i
                jm = j
                km = k
             end if
!!$             data = data + u_next(i,j,k)
          end do                   ! DU COUP LEVEL A ETE AJOUTE POUR L'INSTANT
       enddo
    end do
    if(level(im,jm) > 0) then
       nom = "eau"
    else
       nom = "air"
    end if
    print*, "Maximum de la vitesse atteint en (",im,",",jm,",",km,") = ",u_next(im,jm,km)," dans l'",nom
!!$    print*, data
!!$    read*,

    u = u_next
    p = p_next

    deallocate(laplace_u,u_grad_u)

  end subroutine projection_method




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !méthode de projection de Chorin diphasique
  subroutine projection_method_diphasique(u,p,rho,rho_centre,mu,g,dt,dx,level,A,ipvt)
    implicit none

     ! vitesse au temps n
    real(kind=8), dimension(:,:,:), intent(inout)                                 :: u
    real(kind=8), dimension(:,:), intent(inout)                                   :: p
    real(kind=8), dimension(:,:), intent(in)                                      :: level

    !vitesse au temps n+1
    real(kind=8), dimension(size(u,1),size(u,2),size(u,3))                        :: u_next
    real(kind=8), dimension(size(p,1),size(p,2))                                  :: p_next
    real(kind=8), dimension(size(p,1)*size(p,2))                                  :: p_next_vect

    !vitesse intermédiaire
    real(kind=8), dimension(:,:,:), allocatable                                   :: u_star

    !paramètres physiques
    real(kind=8), dimension(2), intent(in)                                        :: g
    real(kind=8), dimension(:,:), intent(in)                                      :: rho, rho_centre, mu
    real(kind=8), dimension(size(rho_centre,1)*size(rho_centre,2))                :: rho_centre_vect

    !paramètres numériques
    real(kind=8), intent(in)                                                      :: dt, dx 
    integer :: N, i, j

    !matrice de Poisson et second membre
    real(kind=8), dimension(:,:), intent(in)                                      :: A
    integer, dimension(:), intent(in)                                             :: ipvt
    real(kind=8), dimension(:), allocatable                                       :: B

    real(kind=8), dimension(:,:,:), allocatable                                   :: laplace_u , u_grad_u , grad_p

    integer                                                                       :: info, im, jm, km, k

    character(len=3)                                                              :: nom

    real(8)                                                                       :: data, temp, temp2


    N = size(u,1)-1
 
    allocate(laplace_u(2:N,2:N,2),u_grad_u(2:N,2:N,2),grad_p(2:N,2:N,2))
    allocate(u_star(N+1,N+1,2))
    allocate(B(N*N))

    !laplace_u = ( u(3:N+1,2:N,:) - 2*u(2:N,2:N,:) + u(1:N-1,2:N,:) ) / (dx*dx) &
    !     & + ( u(2:N,3:N+1,:) - 2*u(2:N,2:N,:) + u(2:N,1:N-1,:) ) / (dx*dx)
    !ecoulement diphasique ==> terme de viscosité différent
    
    do i = 2 , N
       do j = 2 , N

          !première composante
          temp = 0.25*(u(i,j,2) + u(i+1,j,2) + u(i,j+1,2) + u(i+1,j+1,2)) - &
               & 0.25*(u(i,j,2) + u(i-1,j,2) + u(i,j+1,2) + u(i-1,j+1,2)) + &
               & u(i,j+1,1) - u(i,j,1)

          temp2 = 0.25*(u(i,j,2)+u(i+1,j,2)+u(i,j-1,2)+u(i+1,j-1,2)) - &
               & 0.25*(u(i,j,2)+u(i-1,j-1,2)+u(i-1,j,2)+u(i,j-1,2)) + & 
               & u(i,j,1) - u(i,j-1,1)

          laplace_u(i,j,1) = ((mu(i+1,j) + mu (i,j))*(u(i+1,j,1) - u(i,j,1)) - (mu(i,j) + mu(i-1,j))*(u(i,j,1) - u(i-1,j,1))) &
               & /(rho(i,j)*dx*dx) + &
               & 0.5*((mu(i,j)+mu(i,j+1))*temp - (mu(i,j-1)+mu(i,j))*temp2)/(rho(i,j)*dx*dx)

          !deuxième composante

          temp = u(i+1,j,2) - u(i,j,2) + &
               & 0.25*(u(i,j,1)+u(i+1,j,1)+u(i,j+1,1)+u(i+1,j+1,1)) - &
               & 0.25*(u(i,j,1)+u(i+1,j,1)+u(i,j-1,1)+u(i+1,j-1,1))
          
          temp2 = u(i,j,2) - u(i-1,j,2) + &
               & 0.25*(u(i,j,1)+u(i-1,j,1)+u(i,j+1,1)+u(i-1,j+1,1)) - &
               & 0.25*(u(i,j,1)+u(i-1,j,1)+u(i,j-1,1)+u(i-1,j-1,1))

          laplace_u(i,j,2) = ((mu(i,j+1)+mu(i,j))*(u(i,j+1,2)-u(i,j,2))-(mu(i,j)+mu(i,j-1))*(u(i,j,2)-u(i,j-1,2)))&
               & /(rho(i,j)*dx*dx) + &
               & 0.5*((mu(i+1,j)+mu(i,j))*temp - (mu(i-1,j)+mu(i,j))*temp2)/(rho(i,j)*dx*dx)


       end do
    end do

    !terme d'inertie
    do i = 2 , N
       do j = 2 , N

          u_grad_u(i,j,:) = u(i,j,1)*( u(i+1,j,:) - u(i-1,j,:) ) / (2*dx) &
               & + u(i,j,2)*( u(i,j+1,:) - u(i,j-1,:) ) / (2*dx)         

       end do
    end do

    u_star = 0
    u_star(:,N+1,1) = 0.

    do i = 1, size(u_star,3)
       u_star(2:N,2:N,i) = u(2:N,2:N,i) + dt*( g(i) + laplace_u(:,:,i) - u_grad_u(:,:,i) )
    end do

    !résolution problème de Poisson 
    !ecoulement diphasique ==> probleme différent
    ! fonction mat_vect_diphasique( X, dx, rho) avec rho aux centres des mailles sous forme vecteur taille N*N

    !remplissage second membre
    do i = 1 , N
       do j = 1 , N

          B(i+N*(j-1)) = 1./(2*dt*dx)*(u_star(i+1,j+1,1)+u_star(i+1,j,1)-u_star(i,j+1,1)-u_star(i,j,1) &
               & +u_star(i,j+1,2)+u_star(i+1,j+1,2)-u_star(i+1,j,2)-u_star(i,j,2))

       end do
    end do

    !gradient conjugué
    rho_centre_vect = vect1(rho_centre)
!!$    B = condi_diphasique(B*dx*dx*0.5d0, rho_centre_vect)
    B = B*dx*dx*0.5d0

    P_next_vect = vect1(P_next)
    call grad_conj_diphasique(P_next_vect,B,dx,rho_centre_vect)
    P_next = unvect1(remoy(P_next_vect,Patm))

    do i = 2 , N

       do j = 2 , N

          grad_p(i,j,1) =  p_next(i,j)+p_next(i,j-1)-p_next(i-1,j)-p_next(i-1,j-1)
          grad_p(i,j,2) =  p_next(i,j)+p_next(i-1,j)-p_next(i-1,j-1)-p_next(i,j-1)

       end do

    end do

    u_next = 0.
    u_next(:,N+1,1) = 0.

    do i = 1, size(u_next,3)
       u_next(2:N,2:N,i) = u_star(2:N,2:N,i) - (dt/(2*dx*rho(2:N,2:N)))*grad_p(2:N,2:N,i)
    end do

    im = 1
    jm = 1
    km = 1
    do i = 1, N
       do j = 1, N
          do k = 1, 2
             if(abs(u_next(i,j,k)) > abs(u_next(im,jm,km))) then
                im = i
                jm = j
                km = k
             end if
!!$             data = data + u_next(i,j,k)
          end do                   ! DU COUP LEVEL A ETE AJOUTE POUR L'INSTANT
       enddo
    end do
    if(level(im,jm) > 0.5) then
       nom = "air"
    else
       nom = "eau"
    end if
    print*, "Maximum de la vitesse atteint en (",im,",",jm,",",km,") = ",u_next(im,jm,km)," dans l'",nom
!!$    print*, data
!!$    read*,

    u = u_next
    p = p_next

    deallocate(laplace_u,u_grad_u)


  end subroutine projection_method_diphasique

  function remoy(X, moy) result(Y)
    
    implicit none

    real(8), dimension(:), intent(in) :: X
    real(8), dimension(size(X)) :: Y
    real(8), intent(in) :: moy

    Y = X - sum(X)/size(X) + moy

  end function remoy

end module projection_methodmod
