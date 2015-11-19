program main

  use transportmod
  use modmainmod
  
  implicit none

  !--------------------------------
  integer, parameter :: N = 100   !
  integer, parameter :: tmax = 1  !
  real(8), parameter :: cfl = 0.9 !
  !--------------------------------

  real(8), dimension(N+1,N+1,2) :: noeuds, vitesses
  real(8), dimension(N,N,2) :: centres
  real(8), dimension(N+1,N+1) :: level
  real(8), dimension(N,N) :: pressions
  real(8) :: dx, dt
  integer :: i, j, Niter, ci, di, ui, k
  character(len=1) :: c, d, u
  
  dx = 1./N
  
  call initcoord(noeuds, centres)

  vitesses = 0.3
  pressions = 1.
  
  do i = 1, N+1
     do j = 1, N+1
        level(i,j) = exp(-((float(i)-(N+1)/2)**2+(float(j)-(N+1)/2)**2)/(N+1))
     end do
  end do

  dt = cfl*dx/maxval(vitesses)
  Niter = ceiling(tmax/dt)

  call write(0, noeuds, level)

  do k = 1, Niter

     call transport(noeuds, vitesses, level, dt, dx)

     call write(k, noeuds, level)

  end do

contains

end program main
