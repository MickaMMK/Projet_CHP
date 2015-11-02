program main
  
  implicit none

  type point
     real(8), dimension(2) :: coord
     type(point), pointer :: pointNE, pointNO, pointSE, pointSO
  end type point

  type vecteur
     real(8), dimension(2) :: coord
     type(point), pointer :: point
  end type vecteur

  type scalaire
     real(8) :: valeur
     type(point), pointer :: point
  end type scalaire

  !----------------------------
  integer, parameter :: N=20  !
  !----------------------------

  type(point), dimension(N,N), target :: centres
  type(point), dimension(N+1,N+1), target :: noeuds
  integer :: i, j

  do i = 1, N+1
     do j = 1, N+1
        noeuds(i,j)%coord = (/float(i-1)/N, float(j-1)/N/)
     end do
  end do

  do i = 1, N
     do j = 1, N
        centres(i,j)%coord = (/(float(i-1)+0.5)/N, (float(j-1)+0.5)/N/)
     end do
  end do

  do i = 1, N
     do j = 1, N
        centres(i,j)%pointNO => noeuds(i,j+1)
        centres(i,j)%pointNE => noeuds(i+1,j+1)
        centres(i,j)%pointSE => noeuds(i+1,j)
        centres(i,j)%pointSO => noeuds(i,j)
     end do
  end do

  do i = 2, N
     do j = 2, N
        noeuds(i,j)%pointNO => centres(i-1,j)
        noeuds(i,j)%pointNE => centres(i,j)
        noeuds(i,j)%pointSE => centres(i,j-1)
        noeuds(i,j)%pointSO => centres(i-1,j-1)
     end do
  end do

  open(unit=12, file="noeuds.dat", status="replace")
  open(unit=13, file="centres.dat", status="replace")
  
  do i = 1, N+1
     do j = 1, N+1
        write(12,*) noeuds(i,j)%coord
     end do
  end do

  do i = 1, N
     do j = 1, N
        write(13,*) centres(i,j)%coord
     end do
  end do

  close(12)
  close(13)

end program main
