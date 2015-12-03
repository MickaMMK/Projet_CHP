module transportmod

  use eulermod
  use locatemod
  use interpmod

!!$%%%%%%% COMMENTAIRES %%%%%%%!!
!!$Résolution de l'équation de transport d(f)/dt + v.grad(f)=0
!!$
!!$ ATTENTION, CA DEVRAIT PAS COMPILER: VERIFIER SYNTAXE ET TYPES
!!$
!!$%%%%%%% FIN COMMENTAIRES %%%%%%%!!
  implicit none

contains
  subroutine transport(noeuds, vitesses, level, dt, dx)

    implicit none

    real(8), dimension(:,:,:), intent(in) :: noeuds, vitesses
    real(8), dimension(:,:), intent(inout) :: level
    real(8), intent(in) :: dt, dx
    real(8), dimension(size(level(:,1)),size(level(1,:))) :: old_level
    real(8), dimension(2) :: coord
    real(8), dimension(4,3) :: sommets
    integer, dimension(4,2) :: indices
    integer :: i,j,k,N
    real(8) :: val
    integer :: methode           !Choisir la méthode d'interpolation. 0 = interpolation de base / 1 = autre interpolation 

    methode = 0
    old_level = level
    N = size(level(:,1))-1

    !boucle sur les noeuds
    do i = 1, N+1
       do j = 1, N+1

          !position au temps précedent
          coord = noeuds(i,j,:)
          call euler(-1.*vitesses(i,j,:), coord, dt) !Modifie coord
          if(coord(1) < 1e-6 .OR. coord(2) < 1e-6 .OR. coord(1) > 1-1e-6 .OR. coord(2) > 1-1e-6) then
          else
             call locate(coord, dx, indices) !Récupère les indices des sommets
             !valeur au temps précédent
             do k = 1, 4
                sommets(k,1:2) = noeuds(indices(k,1), indices(k,2),:)  !sommets contient pour les 4 sommets
                sommets(k,3) = old_level(indices(k,1), indices(k,2)) !les coord x, y et la valeur du level
             end do
             !call interp(sommets, coord, val, methode) !Interpolation de la valeur du level
             call noyau_interp(sommets,coord,dx,val)
             level(i,j) = val
          end if

       end do
    end do

  end subroutine transport

end module transportmod
