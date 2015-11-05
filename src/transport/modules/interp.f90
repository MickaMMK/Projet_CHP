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

  function poids(a, b, x, y) result(p)
    
!!$%%%%%%% COMMENTAIRES %%%%%%%!!
!!$
!!$ Pondération : p(a,b) = 1/sqrt((a-x)**2+(b-y)**2), où (x,y) sont les coordonnées du point où l'on souhaite connaitre
!!$
!!$ la valeur de phi
!!$
!!$%%%%%%% FIN COMMENTAIRES %%%%%%%!!
    
    implicit none
    
    real(8), intent(in) :: a,b,x,y
    real(8) :: p
    
    p = 1.0/sqrt((a-x)**2+(b-y)**2)

  end function poids
  


  subroutine interp(sommets, coord, valeur)

!!$%%%%%%% COMMENTAIRES %%%%%%%!!
!!$
!!$ Interpolation bilinéaire dans un carré avec pondération donnée par la distance du point à chaque sommet
!!$
!!$ La fonction poids calcule cette pondération
!!$
!!$%%%%%%% FIN COMMENTAIRES %%%%%%%!!

    implicit none
    
    !Déclaration des variables
    
    real(8),dimension(4,3),intent(in) :: sommets   ! Tableau 4*3. Chaque case contient les informations sur le sommet(x_pos, y_pos,valeur de phi à ce point)
    real(8),dimension(2),intent(in) :: coord       ! Coordonnées du point où l'on veut savoir la valeur de phi
    real(8),intent(out) :: valeur
    
    integer :: sommet
    real(8) :: ponderation,p

    valeur = 0
    ponderation = 0

    do sommet = 1, 4

       p= poids(sommets(sommet,1),sommets(sommet,2),coord(1),coord(2))

       valeur = valeur + p*sommets(sommet,3)

       ponderation = ponderation + p

    end do
    
    valeur = valeur/ponderation



  end subroutine interp


end module interpmod
