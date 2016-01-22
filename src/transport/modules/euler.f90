module eulermod
!!$%%%%%%% COMMENTAIRES %%%%%%%!!
!!$ Méthode d'Euler pour situer au temps n-1 la valeur sur un noeud au temps n
!!$ de la fonction cherchée 
!!$%%%%%%% FIN COMMENTAIRES %%%%%%%!!
  implicit none

   !---------------------------------------------------------------------------
   !> @brief
   !> 
   !
   !> @author 
   !> Corentin PRIGENT, engineer student at ENSEIRB-MATMECA, Bordeaux, FR.
   !
   !> @param[in] dx, N      
   !> @param[out] A      
   !> @return The Poisson's matrix
   !--------------------------------------------------------------------------- 
 
contains
  subroutine euler(vitesse, coord, dt)
    implicit none
    
    !Déclaration des variables
    real(8),dimension(2),intent(in) :: vitesse
    real(8),dimension(2),intent(inout) :: coord
    real(8),intent(in) :: dt

    coord = coord + dt*vitesse

  end subroutine euler


end module eulermod
