module transportmod
!!$%%%%%%% COMMENTAIRES %%%%%%%!!
!!$Résolution de l'équation de transport d(f)/dt + v.grad(f)=0
!!$
!!$ ATTENTION, CA DEVRAIT PAS COMPILER: VERIFIER SYNTAXE ET TYPES
!!$
!!$%%%%%%% FIN COMMENTAIRES %%%%%%%!!
  implicit none

contains
  subroutine transport(noeuds, vitesses, level, dt)
    real(8), dimension(N+1,N+1,2), intent(in) :: noeuds, vitesses
    real(8), dimension(N+1,N+1), intent(inout) :: level
    real(8), intent(in) :: dt
    real(8), dimension(N+1,N+1) :: old_level
    real(8), dimension(2) :: coord
    real(8), dimension(4,3) :: sommets
    integer, dimension(4,2) :: indices
    integer::i,j,k

    old_level = level

    !boucle sur les noeuds
    do i = 1, N+1
       do j = 1, N+1

          !position au temps précedent
          coord = noeuds(i,j)
          call euler(-1.*vitesses(i,j,:), coord, dt) !Modifie coord
          call find_noeuds(coord, indices) !Récupère les indices des sommets
          !valeur au temps précédent
          do k = 1, 4
             sommets(k,1:2) = noeuds(indices(k,1), indices(k,2))  !sommets contient pour les 4 sommets
             sommets(k,3) = old_level(indices(k,1), indices(k,2)) !les coord x, y et la valeur du level
          end do
          call interp(sommets, coord, val) !Interpolation de la valeur du level
          level(i,j) = val

       end do
    end do

  end subroutine transport

end module transportmod
