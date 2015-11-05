program main

  use transportmod
  
  implicit none

  !--------------------------------
  integer, parameter :: N = 5     !
  integer, parameter :: tmax = 1  !
  real(8), parameter :: cfl = 0.9 !
  !--------------------------------

  real(8), dimension(N+1,N+1,2) :: noeuds, vitesses
  real(8), dimension(N,N,2) :: centres
  real(8), dimension(N+1,N+1) :: level
  real(8), dimension(N,N) :: pressions
  real(8) :: dx, dt
  integer :: i, j, Niter
  
  dx = 1./N
  do i = 1, N+1
     do j = 1, N+1
        noeuds(i,j,:) = (/(i-1)*dx, (j-1)*dx/)
        vitesses(i,j,:) = (/0.4,0.4/)
        if(sqrt((noeuds(i,j,1)-0.5)**2+(noeuds(i,j,2)-0.5)**2) < 0.2) then
           level(i,j) = 1.
        else
           level(i,j) = 0.
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

  open(unit=22, file="noeuds.dat", status="replace")
  open(unit=23, file="centres.dat", status="replace")
  
  do i = 1, N+1
     do j = 1, N+1
        write(22,*) noeuds(i,j,:), level(i,j)
     end do
     write(22,*)
  end do

  do i = 1, N
     do j = 1, N
        write(23,*) centres(i,j,:)
     end do
     write(23,*)
  end do

  close(22)
  close(23)

  do i = 1, Niter
     call transport(noeuds, vitesses, level, dt, dx)
  end do

  open(unit=12, file="noeudsf.dat", status="replace")
  open(unit=13, file="centresf.dat", status="replace")
  
  do i = 1, N+1
     do j = 1, N+1
        write(12,*) noeuds(i,j,:), level(i,j)
     end do
     write(12,*)
  end do

  do i = 1, N
     do j = 1, N
        write(13,*) centres(i,j,:)
     end do
     write(13,*)
  end do

  close(12)
  close(13)



contains

end program main
