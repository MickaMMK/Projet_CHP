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
  


  subroutine interp(sommets, coord, valeur, methode)

!!$%%%%%%% COMMENTAIRES %%%%%%%!!
!!$
!!$ Methode = 0 : Interpolation bilinéaire dans un carré avec pondération donnée par la distance du point à chaque sommet.
!!$
!!$               La fonction poids calcule cette pondération.
!!$
!!$ Methode 1 : Interpolation de la forme phi(x,y) = aX+bY+cXY+d,
!!$
!!$             avec X=x-x4 et Y=y-y4 au point de coordonnées (x,y).
!!$
!!$%%%%%%% FIN COMMENTAIRES %%%%%%%!!

    implicit none
    
    !Déclaration des variables
    
    real(8),dimension(4,3),intent(in) :: sommets   ! Tableau 4*3. Chaque case contient les informations sur le sommet(x_pos, y_pos,valeur de phi à ce point)
    real(8),dimension(2),intent(in) :: coord       ! Coordonnées du point où l'on veut savoir la valeur de phi
    integer,intent(in) :: methode
    real(8),intent(out) :: valeur
    
    integer :: sommet
    real(8) :: ponderation,p

    real(8) :: dx,dy
    real(8) :: a,b,c,d
    
    real(8) :: X,Y

    if (methode ==0) then


       valeur = 0
       ponderation = 0

       do sommet = 1, 4

          p= poids(sommets(sommet,1),sommets(sommet,2),coord(1),coord(2))

          valeur = valeur + p*sommets(sommet,3)

          ponderation = ponderation + p

       end do

       valeur = valeur/ponderation


    elseif (methode == 1) then


       ! Définition de dx et dy, les dimensions du carré

       dx=abs(sommets(3,1)-sommets(4,1)) !x3-x4
       dy=abs(sommets(1,2)-sommets(4,2)) !y1-y4


       !Définition des a,b,c et d

       a=(sommets(3,3)-sommets(4,3))/dx
       b=(sommets(1,3)-sommets(4,3))/dy
       c=((sommets(2,3)+sommets(4,3))-sommets(1,3)-sommets(3,3))/(dx*dy)
       d=sommets(4,3)

       !Valeur de phi

       X=coord(1)-sommets(4,1)
       Y=coord(2)-sommets(4,2)

       valeur = a*X+b*Y+c*X*Y+d


    end if

  end subroutine interp
  


end module interpmod
