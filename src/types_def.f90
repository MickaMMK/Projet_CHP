program main
  
  implicit none

  !----------------------------
  integer, parameter :: N=100 !
  !----------------------------

  real(8), dimension(N+1,N+1,2) :: coord, vitesses
  real(8), dimension(N,N,2) :: centres
  real(8), dimension(N+1,N+1) :: level
  real(8), dimension(N,N) :: pressions
  real(8) :: dx
  integer :: i, j

  dx = 1/N

  do i = 1, N+1
     do j = 1, N+1
        coord(i,j,:) = (/(i-1)*dx, (j-1)*dx/)
        vitesses(i,j,:) = (/1,1/)
        if(sqrt((coord(i,j,1)-0.5)**2+(coord(i,j,2)-0.5)**2) < 0.2) then
           level(i,j) = 1.
        else
           level(i,j) = 0.
        endif
     end do
  end do

  do i = 1, N
     do j = 1, N
        centres(i,j,:) = (/((i-1)+0.5)*dx, (j-1)+0.5)*dx/)
        pressions(i,j) = 1.
     end do
  end do


  open(unit=12, file="noeuds.dat", status="replace")
  open(unit=13, file="centres.dat", status="replace")
  
  do i = 1, N+1
     do j = 1, N+1
        write(12,*) coord(i,j,:), level(i,j)
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
