module transportmod

  use eulermod
  use locatemod
  use interpmod
  use vectmod

  implicit none

contains
  subroutine transport_level(noeuds, vitesses, level, dt, dx)

    implicit none

    real(8), dimension(:,:,:), intent(in) :: noeuds, vitesses
    real(8), dimension(:,:), intent(inout) :: level
    real(8), intent(in) :: dt, dx
    real(8), dimension(size(noeuds,1)*size(noeuds,2),size(noeuds,3)) :: vect_noeuds
    real(8), dimension(size(level,1)*size(level,2)) :: old_level
    real(8), dimension(2) :: coord
    integer :: i,j,k,N
    real(8) :: val

    N = size(level,1)-1
    old_level = vect1(level)
    vect_noeuds = vect2(noeuds)

    !boucle sur les noeuds
    do i = 1, N+1
       do j = 1, N+1

          !position au temps précedent
          coord = noeuds(i,j,:)
          call euler(-1.*vitesses(i,j,:), coord, dt) !Modifie coord
          if(coord(1) < 0 .OR. coord(2) < 0 .OR. coord(1) > 1 .OR. coord(2) > 1) then
          else
             call noyau_interp(vect_noeuds,old_level,coord,dx,val)
             level(i,j) = val
          end if

       end do
    end do

  end subroutine transport_level

  subroutine transport_particules(particules, vitesses_particules, dt, dx)

    implicit none

    real(8), dimension(:,:), intent(inout) :: particules
    real(8), dimension(:,:), intent(in) :: vitesses_particules
    real(8), intent(in) :: dt, dx
    integer :: i,j,Nm

    Nm = size(particules,1)
    do i = 1, Nm
          call euler(vitesses_particules(i,:), particules(i,:), dt)
    end do
    

  end subroutine transport_particules

  subroutine transport_level_EL(noeuds, vitesses, level, dt, dx, coordlagr, levellagr)

    implicit none

    real(8), dimension(:,:,:), intent(in) :: noeuds, vitesses
    real(8), dimension(:,:), intent(in) :: coordlagr
    real(8), dimension(:,:), intent(inout) :: level
    real(8), dimension(:), intent(in) :: levellagr
    real(8), intent(in) :: dt, dx
    real(8), dimension(2) :: coord
    real(8), dimension(size(noeuds,1)*size(noeuds,2)+size(levellagr),2) :: coordtot
    real(8), dimension(size(noeuds,1)*size(noeuds,2)+size(levellagr)) :: leveltot
    integer :: i,j,k,N
    real(8) :: val

    coordtot(1:size(noeuds,1)*size(noeuds,2),:) = vect2(noeuds)
    coordtot(size(noeuds,1)*size(noeuds,2)+1:size(coordtot),:) = coordlagr

    leveltot(1:size(noeuds,1)*size(noeuds,2)) = vect1(level)
    leveltot(size(noeuds,1)*size(noeuds,2)+1:size(coordtot)) = levellagr

    !boucle sur les noeuds
    do i = 1, size(noeuds,1)
       do j = 1, size(noeuds,2)

          !position au temps précedent
          coord = noeuds(i,j,:)
          call euler(-1.*vitesses(i,j,:), coord, dt) !Modifie coord
          if(coord(1) < 1e-6 .OR. coord(2) < 1e-6 .OR. coord(1) > 1-1e-6 .OR. coord(2) > 1-1e-6) then
          else
             call noyau_interp(coordtot,leveltot,coord,dx,val)
             level(i,j) = val
          end if

       end do
    end do

  end subroutine transport_level_EL

end module transportmod
