program main

  use transportmod
  use modmainmod
  use projection_methodmod
  use remplissage_poissonmod
  
  implicit none

  !--------------------------------
  integer, parameter :: N = 50    !
  integer, parameter :: tmax = 5  !
  real(8), parameter :: cfl = 0.9 !
  !--------------------------------

  real(8), parameter :: pi = 3.1415926535897932384626433842795028841971693993

  real(8), dimension(N+1,N+1,2) :: noeuds, vitesses
  real(8), dimension(:,:), allocatable :: particules, vitesses_particules
  real(8), dimension(N,N,2) :: centres
  real(8), dimension(N+1,N+1) :: level, rho, nu
  real(8), dimension(:), allocatable :: level_particules
  real(8), dimension(N,N) :: pressions, rho_centre
  real(8), dimension(N*N,N*N) :: A
  integer, dimension(N*N) :: ipvt
  integer :: info
  real(8), dimension(2) :: g
  real(8) :: dx, dt, rho_air, rho_eau, nu_air, nu_eau
  integer :: i, j, Niter, ci, di, ui, k, meth, nbp
  character(len=1) :: c, d, u
  
  dx = 1./N

  rho_air = 1000. !1.
  rho_eau = 1000. !1000.
  nu_air = 0.9E-2 !15-6
  nu_eau = 0.9E-2 !0.9E-6
  g = (/0.,-9.81/)

  print*, "Choix de la méthode :"
  print*, "1 - Eulerien"
  print*, "2 - Lagrangien"
  print*, "3 - Eulerien-Lagrangien"
  read*, meth
  if(meth < 1 .OR. meth > 3) then
     print*, "Erreur dans le choix de la méthode"
     print*, "==========="
     print*, "== ABORT =="
     print*, "==========="
     stop
  else if(meth == 2) then
     nbp = (N-1)*(N-1)
  else if(meth == 3) then
     nbp = 2000
  end if

  allocate(particules(nbp,2), vitesses_particules(nbp,2), level_particules(nbp))
  
  call initcoord(noeuds, centres)

  vitesses = 0.
!!$  vitesses(2:N,2:N,:) = 0.3

!!$  do i = 2, N
!!$     do j = 2, N
!!$        vitesses(i,j,:) = (/-(noeuds(i,j,1)-0.5)**2+0.25,-(noeuds(i,j,2)-0.5)**2+0.25/)
!!$        print*, vitesses(i,j,:)
!!$     end do
!!$  end do
  pressions = 1.013D5
  
  do i = 1, N+1
     do j = 1, N+1
        if(sqrt((noeuds(i,j,1)-0.5)**2+(noeuds(i,j,2)-0.5)**2) < 0.2) then
           level(i,j) = 0.
           vitesses(i,j,2) = -2.
        else
           level(i,j) = 1.
        end if
     end do
  end do

  if(meth == 2) then
     particules = vect2(noeuds(2:N,2:N,:))
     level_particules = vect1(level(2:N,2:N))
  else if(meth == 3) then
     do i = 1, nbp/2
        particules(i,:) = (/0.5+0.19*cos(4*pi*i/nbp), 0.5+0.19*sin(4*pi*i/nbp)/)
        level_particules(i) = 0.
        particules(nbp/2+i,:) = (/0.5+0.21*cos(4*pi*i/nbp), 0.5+0.21*sin(4*pi*i/nbp)/)
        level_particules(nbp/2+i) = 1.
     end do
  end if

  print*, maxval(abs(vitesses))
  dt = cfl*min(dx/maxval(abs(vitesses)),dx*dx/(4*max(nu_air,nu_eau)),2*min(nu_air,nu_eau)/maxval(abs(vitesses))**2)
  print*, "dt = ",dt
  Niter = ceiling(tmax/dt)

  call write("level", 0, vect2(noeuds), vect1(level))
  call write("vitesses_x", 0, vect2(noeuds), vect1(vitesses(:,:,1)))
  call write("vitesses_y", 0, vect2(noeuds), vect1(vitesses(:,:,2)))
  call write("pressions", 0, vect2(centres), vect1(pressions))

  call remplissage_poisson(A,dx,N)
  call DGETRF(N*N, N*N, A, N*N, ipvt, info)

  do k = 1, Niter

     print*, "Itération ",k," sur ",Niter
  
     rho = int(level+0.5)*rho_air + (1-int(level+0.5))*rho_eau
     nu = int(level+0.5)*nu_air + (1-int(level+0.5))*nu_eau
     rho_centre = (rho(1:N,1:N)+rho(1:N,2:N+1)+rho(2:N+1,1:N)+rho(2:N+1,2:N+1))/4

     if(meth == 1) then

        !==============================================================================================
        !========================================== EULERIEN ==========================================
        
        print*, 'Projection'
        call projection_method(vitesses, pressions, rho, rho_centre, nu, g, dt, dx, level, A, ipvt)

        print*, 'Transport'
        call transport_level(noeuds, vitesses, level, dt, dx)

        !print*, "somme des lvl = ",sum(1-level)
        !print*, "somme des lvl après valeur fixée à 0 ou 1 = ",sum(int(1-level+0.5))

        print*, 'Ecriture'
        call write("level", k, vect2(noeuds), vect1(level))
        call write("vitesses_x", k, vect2(noeuds), vect1(vitesses(:,:,1)))
        call write("vitesses_y", k, vect2(noeuds), vect1(vitesses(:,:,2)))
        call write("pressions", k, vect2(centres), vect1(pressions))

        !dt = cfl*min(dx/maxval(vitesses),dx*dx/(2*max(nu_air,nu_eau)))

        !==============================================================================================

     else if(meth == 2) then
        
        !==============================================================================================
        !========================================= LAGRANGIEN =========================================
        
        print*, 'Projection'
        call projection_method(vitesses, pressions, rho, rho_centre, nu, g, dt, dx, level, A, ipvt)

        do i = 1, (N-1)*(N-1)
           call noyau_interp(vect2(noeuds), vect1(vitesses(:,:,1)), particules(i,:), dx, vitesses_particules(i,1))
           call noyau_interp(vect2(noeuds), vect1(vitesses(:,:,2)), particules(i,:), dx, vitesses_particules(i,2))
        end do

        call transport_particules(particules, vitesses_particules, dt, dx)

        do i = 2, N
           do j = 2, N
              call noyau_interp(particules, level_particules, noeuds(i,j,:), dx, level(i,j))
           end do
        end do

        if(mod(k,50) == 0) then
           particules = vect2(noeuds(2:N,2:N,:))
           vitesses_particules = vect2(vitesses(2:N,2:N,:))
           level_particules = vect1(level(2:N,2:N))
           !level_particules = int(level(2:N,2:N)+0.5)
        end if

        call write("level_particules", k, particules, level_particules)
!!$     call write("vitesses_x_particules", k, particules, vitesses_particules(:,1))
!!$     call write("vitesses_y_particules", k, particules, vitesses_particules(:,2))
!!$
!!$     call write("level", k, vect2(noeuds), vect1(level))
!!$     call write("vitesses_x", k, vect2(noeuds), vect1(vitesses(:,:,1)))
!!$     call write("vitesses_y", k, vect2(noeuds), vect1(vitesses(:,:,2)))
!!$     call write("pressions", k, vect2(centres), vect1(pressions))

        !==============================================================================================

     else if(meth == 3) then
        
        !==============================================================================================
        !===================================== EULERIEN-LAGRANGIEN ====================================
        
        print*, 'Projection'
        call projection_method(vitesses, pressions, rho, rho_centre, nu, g, dt, dx, level, A, ipvt)

        do i = 1, nbp
           call noyau_interp(vect2(noeuds), vect1(vitesses(:,:,1)), particules(i,:), dx, vitesses_particules(i,1))
           call noyau_interp(vect2(noeuds), vect1(vitesses(:,:,2)), particules(i,:), dx, vitesses_particules(i,2))
        end do

        print*, 'Transport'
        call transport_level_EL(noeuds, vitesses, level, dt, dx, particules, level_particules)
        call transport_particules(particules, vitesses_particules, dt, dx)

        call write("level_particules", k, particules, level_particules)
        call write("vitesses_x_particules", k, particules, vitesses_particules(:,1))
        call write("vitesses_y_particules", k, particules, vitesses_particules(:,2))

        call write("level", k, vect2(noeuds), vect1(level))
        call write("vitesses_x", k, vect2(noeuds), vect1(vitesses(:,:,1)))
        call write("vitesses_y", k, vect2(noeuds), vect1(vitesses(:,:,2)))
        call write("pressions", k, vect2(centres), vect1(pressions))

        !==============================================================================================

     else
        
        print*, "Putain mais comment t'es arrivé là avec une mauvaise méthode ?!"

     end if

  end do

end program main
