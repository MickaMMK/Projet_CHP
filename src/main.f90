program main

  use transportmod
  use modmainmod
  use navierstokes
  
  implicit none

  !--------------------------------
  integer, parameter :: N = 50     !
  integer, parameter :: tmax = 1  !
  real(8), parameter :: cfl = 0.9 !
  !--------------------------------

  real(8), dimension(N+1,N+1,2) :: noeuds, vitesses
  real(8), dimension(N,N,2) :: centres
  real(8), dimension(N+1,N+1) :: level, rho, nu
  real(8), dimension(N,N) :: pressions, rho_centre
  integer, dimension(N) :: ipvt
  integer :: info
  real(8), dimension(2) :: g
  real(8) :: dx, dt, rho_air, rho_eau, nu_air, nu_eau
  integer :: i, j, Niter, ci, di, ui, k
  character(len=1) :: c, d, u
  
  dx = 1./N

  rho_air = 1.
  rho_eau = 1000.
  nu_air = 15E-6
  nu_eau = 0.9E-6
  g = (/0.,-9.81/)
  
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
        else
           level(i,j) = 1.
        end if
     end do
  end do

  dt = cfl*dx/10!maxval(vitesses)
  print*, "dt = ",dt
  Niter = ceiling(tmax/dt)

  call write(0, noeuds, level)

  do k = 1, Niter
  
     rho = int(level+0.5)*rho_air + (1-int(level+0.5))*rho_eau
     nu = int(level+0.5)*nu_air + (1-int(level+0.5))*nu_eau
     rho_centre = (rho(1:N,1:N)+rho(1:N,2:N+1)+rho(2:N+1,1:N)+rho(2:N+1,2:N+1))/4
     !print*, pressions

     print*, 'Projection'
     call projection_method(vitesses, pressions, rho, rho_centre, nu, g, dt, dx, level)

     print*, 'Transport'
     call transport(noeuds, vitesses, level, dt, dx)

     print*, 'Ecriture'
     call write(k, noeuds, rho)

  end do

contains

end program main
