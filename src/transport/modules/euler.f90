module eulermod
!!$%%%%%%% COMMENTAIRES %%%%%%%!!
!!$ Méthode d'Euler pour situer au temps n-1 la valeur sur un noeud au temps n
!!$ de la fonction cherchée 
!!$
!!$ ATTENTION, CA DEVRAIT PAS COMPILER: VERIFIER SYNTAXE ET TYPES
!!$
!!$%%%%%%% FIN COMMENTAIRES %%%%%%%!!
  implicit none

contains
  subroutine euler(vitesse, coord, dt)
    implicit none
    
    !Déclaration des variables
    real(8),dimension(2),intent(in) :: vitesse
    real(8),dimension(2),intent(inout) : coord
    real(8),intent(in) :: dt

    coord = coord + dt*vitesse

  end subroutine euler


end module eulermod
