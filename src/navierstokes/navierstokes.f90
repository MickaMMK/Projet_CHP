module navierstokes
  use grad_conj
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

    A(N*N,N*N) = dx*dx

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

  end subroutine remplissage_poisson

  !méthode de projection de Chorin
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
    integer, dimension(:), intent(in)          :: ipvt
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

    B(N*N) = 1.013D5

    do i = 1, N
       do j = 1, N
          p_next_vect(i+(j-1)*N) = p_next(i,j)
       end do
    end do


    !gradient conjugué

    B = B*dx*dx

!!$    call grad_conj_opt(P_next_vect,B,dx)

!!$    call DGETRF(N*N, N*N, A, N*N, ipvt, info)
    call DGETRS('N', N*N, 1, A, N*N, ipvt, B, N*N, info)
    P_next_vect = B    !calcul vitesse au temps n+1

    do i = 1, N
       do j = 1, N
           p_next(i,j) = p_next_vect(i+(j-1)*N)
       end do
    end do

    do i = 2 , N

       do j = 2 , N

          grad_p(i,j,1) =  p_next(i,j)+p_next(i,j-1)-p_next(i-1,j)-p_next(i-1,j-1)
          grad_p(i,j,2) =  p_next(i,j)+p_next(i-1,j)-p_next(i-1,j-1)-p_next(i,j-1)

       end do

    end do

    u_next = 0

    do i = 1, size(u_next,3)
       u_next(2:N,2:N,i) = u_star(2:N,2:N,i) - (dt/(2*dx*rho(2:N,2:N)))*grad_p(2:N,2:N,i)
    end do
   
    im = 1
    jm = 1
    km = 1
    do i = 1, N
       do j = 1, N
          do k = 1, 2
!!$             if(abs(u_next(i,j,k)) > abs(u_next(im,jm,km))) then
!!$                im = i
!!$                jm = j
!!$                km = k
!!$             end if
             data = data + u_next(i,j,k)
          end do                   ! DU COUP LEVEL A ETE AJOUTE POUR L'INSTANT
       enddo
    end do
!!$    if(level(im,jm) > 0.5) then
!!$       nom = "air"
!!$    else
!!$       nom = "eau"
!!$    end if
!!$    print*, "Maximum de la vitesse atteint en (",im,",",jm,",",km,") = ",u_next(im,jm,km)," dans l'",nom
    print*, data
!!$    read*,

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
