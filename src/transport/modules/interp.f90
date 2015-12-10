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
    
    real(8),dimension(:,:),intent(in) :: sommets   ! Tableau ?*3. Chaque case contient les informations sur le sommet(x_pos, y_pos,valeur de phi à ce point)
    real(8),dimension(2),intent(in) :: coord       ! Coordonnées du point où l'on veut savoir la valeur de phi
    integer,intent(in) :: methode
    real(8),intent(out) :: valeur
    
    integer :: sommet
    real(8) :: ponderation,p

    real(8) :: dx,dy
    real(8) :: a,b,c,d
    
    real(8) :: X,Y

    if (methode ==0) then


       valeur = 0.
       ponderation = 0.

       do sommet = 1, size(sommets,1)

          p= poids(sommets(sommet,1),sommets(sommet,2),coord(1),coord(2))

          valeur = valeur + p*sommets(sommet,3)

          ponderation = ponderation + p

       end do

       valeur = valeur/ponderation


    elseif (methode == 1 .AND. size(sommets,1) == 4) then


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

    else

       print*, "================================================="
       print*, "Méthode inconnue ou appel de méthode 1 hors carré"
       print*, "================================================="

    end if

  end subroutine interp


  subroutine noyau_interp(noeuds, level, coord, dx, valeur)

!!$%%%%%%% COMMENTAIRES %%%%%%%!!
!!$
!!$ Noyau d'interpolation M'_{4} . CF Inviscid axisymmetrization of an elliptical vortex
!!$
!!$%%%%%%% FIN COMMENTAIRES %%%%%%%!!

    implicit none
    
    real(8),dimension(:,:),intent(in) :: noeuds    ! Tableau ?*2. Chaque case contient les coordonnées
    real(8),dimension(:),intent(in)   :: level     ! Tableau ?. Chaque case contient le level
    real(8),dimension(2),intent(in)   :: coord     ! Coordonnées du point où l'on veut savoir la valeur de phi
    real(8),intent(in)                :: dx        ! Pas d'espace
    real(8),intent(out)               :: valeur

    integer                           :: i, j
    real(8)                           :: ponderation
    real(8)                           :: dist


    ponderation = 0.
    valeur = 0.

    
    do i = 1, size(noeuds,1)
       
       dist = sqrt( (noeuds(i,1)-coord(1))**2 + (noeuds(i,2)-coord(2))**2 )/dx

       if ( dist <= 2.0 ) then

          if ( dist <= 1.0 ) then

             valeur = valeur + level(i)*((2-dist)**3- (4*(1-dist)**3))/6
             !valeur = valeur + level(i,j)*(1.0 - (5.0*dist**2)/2.0 + (3.0*dist**3)/2.0)

             ponderation = ponderation +((2-dist)**3- (4*(1-dist)**3))/6
             !ponderation = ponderation +(1.0 - (5.0*dist**2)/2.0 + (3.0*dist**3)/2.0)

          else

             valeur = valeur + level(i)*((2-dist)**3/6)
             !valeur = valeur + level(i,j)*0.5*(1.0-dist)*(2.0-dist)**2

             ponderation = ponderation + ((2-dist)**3/6)
             !ponderation = ponderation + 0.5*(1.0-dist)*(2.0-dist)**2

          end if

       end if

    end do

    if (ponderation < 1d-8) then  
       valeur = valeur/ponderation
    end if

  end subroutine noyau_interp

  function vect1(array) result(vector)

    implicit none
    
    real(8), dimension(:,:), intent(in) :: array
    real(8), dimension(size(array,1)*size(array,2)) :: vector
    integer :: i, j

    do i = 1, size(array,1)
       do j = 1, size(array,2)
          vector(i+size(array,1)*(j-1)) = array(i,j)
       end do
    end do

  end function vect1

  function vect2(array) result(vector)

    implicit none
    
    real(8), dimension(:,:,:), intent(in) :: array
    real(8), dimension(size(array,1)*size(array,2),size(array,3)) :: vector
    integer :: i, j

    do i = 1, size(array,1)
       do j = 1, size(array,2)
          vector(i+size(array,1)*(j-1),:) = array(i,j,:)
       end do
    end do

  end function vect2



end module interpmod
