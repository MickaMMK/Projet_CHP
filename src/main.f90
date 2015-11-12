program main

  use transportmod
  
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
  do i = 1, N+1
     do j = 1, N+1
        noeuds(i,j,:) = (/(i-1)*dx, (j-1)*dx/)
        vitesses(i,j,:) = (/0.3,0.3/)
        if(sqrt((noeuds(i,j,1)-0.5)**2+(noeuds(i,j,2)-0.5)**2) < 0.2) then
           level(i,j) = 0.
        else
           level(i,j) = 1.
        endif
     end do
  end do

  do i = 1, N
     do j = 1, N
        centres(i,j,:) = (/((i-1)+0.5)*dx, ((j-1)+0.5)*dx/)
        pressions(i,j) = 1.
     end do
  end do

  dt = cfl*dx/maxval(vitesses)
  Niter = ceiling(tmax/dt)

  open(unit=22, file="noeuds000.txt", status="replace")

  write(22,*) '"x", "y", "level"'
  
  do i = 1, N+1
     do j = 1, N+1
        write(22,*) noeuds(i,j,1),", ", noeuds(i,j,2),", ", level(i,j)
     end do
  end do

  close(22)

  do k = 1, Niter

     call transport(noeuds, vitesses, level, dt, dx)

     ci = int(k/100)
     di = modulo(int(k/10),10)
     ui = modulo(k,10)

     write(c,'(I1)') ci
     write(d,'(I1)') di
     write(u,'(I1)') ui

     open(unit=12, file="noeuds"//c//d//u//".txt", status="replace")

     write(12,*) '"x", "y", "level"'

     do i = 1, N+1
        do j = 1, N+1
           write(12,*) noeuds(i,j,1),", ", noeuds(i,j,2),", ", level(i,j)
        end do
     end do
     
     close(12)

  end do

contains

end program main
