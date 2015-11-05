module interpmod
!!$%%%%%%% COMMENTAIRES %%%%%%%!!
!!$ 
!!$ Interpolation dans un carré. On a la valeur d'une fonction phi à chaque sommet du
!!$
!!$ carré et on interpole la valeur de cette fonction sur un point à l'intérieur du carré
!!$ 
!!$
!!$%%%%%%% FIN COMMENTAIRES %%%%%%%!!
  implicit none

contains
  subroutine interp(sommets, coord, valeur)
    implicit none
    
    !Déclaration des variables
    
    real(8),dimension(4,3),intent(in) :: sommets   ! Tableau 4*3. Chaque case contient les informations sur le sommet(x_pos, y_pos,valeur de phi à ce point)
    real(8),dimension(2),intent(in) :: coord       ! Coordonnées du point où l'on veut savoir la valeur de phi
    real(8),intent(out) :: valeur

  end subroutine interp


end module interpmod
