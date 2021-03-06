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
!!$ Noyau d'interpolation M'_{4} . CF Inviscid axisymetrization of an elliptical vortex
!!$
!!$%%%%%%% FIN COMMENTAIRES %%%%%%%!!

    implicit none
    
    real(8),dimension(:,:),intent(in) :: noeuds    ! Tableau ?*2. Chaque case contient les coordonnées
    real(8),dimension(:),intent(in)   :: level     ! Tableau ?. Chaque case contient le level
    real(8),dimension(2),intent(in)   :: coord     ! Coordonnées du point où l'on veut savoir la valeur de phi
    real(8),intent(in)                :: dx        ! Pas d'espace
    real(8),intent(out)               :: valeur

    integer                           :: i, j, compteur
    real(8)                           :: ponderation
    real(8)                           :: dist, temp
    real(8)                           :: u, v      ! u=x/dx; v=y/dy
    


    ponderation = 0.
    valeur = 0.
    compteur = 0

    ! Kernels multipliés
    do i = 1, size(noeuds,1)

       u = abs(noeuds(i,1)-coord(1))/dx
       v = abs(noeuds(i,2)-coord(2))/dx

!!$       !Ligne du haut
!!$       !if ( ( v >= -2.0 ) .and. ( v <= -1.0 ) ) then
!!$       if (abs(v+1.5)<=0.5) then
!!$
!!$          !Carré gauche
!!$          !if ( ( u >= -2.0) .and. ( u <= -1.0) ) then
!!$          if (abs(u+1.5)<=0.5) then
!!$
!!$             temp = ((2-v)**3/6)*((2-u)**3/6)
!!$
!!$             valeur = valeur + level(i)*temp
!!$
!!$             ponderation = ponderation + temp
!!$
!!$          !Carré milieu 
!!$          !elseif ( ( u >= -1.0) .and. ( u <= 1.0) ) then
!!$          elseif (abs(u)<1) then
!!$
!!$             temp = ((2-v)**3/6)*((2-u)**3/6 - 4*(1-u)**3/6)
!!$
!!$             valeur = valeur + level(i)*temp
!!$
!!$             ponderation = ponderation + temp
!!$
!!$          !Carré droit            
!!$          !elseif ( ( u >= 1.0) .and. ( u <= 2.0) ) then
!!$          elseif (abs(u-1.5)<=0.5) then
!!$
!!$             temp = ((2-v)**3/6)*((2-u)**3/6)
!!$
!!$             valeur = valeur + level(i)*temp
!!$
!!$             ponderation = ponderation + temp
!!$
!!$          end if


       !Ligne du milieu
       !elseif ( ( v >= -1.0 ) .and. ( v <= 1.0 ) ) then
       if (v<=1.) then

!!$          !Carré gauche
!!$          !if ( ( u >= -2.0) .and. ( u <= -1.0) ) then
!!$          if (abs(u+1.5)<=0.5) then
!!$
!!$             temp = ((2-v)**3/6 - 4*(1-v)**3/6)*((2-u)**3/6)
!!$
!!$             valeur = valeur + level(i)*temp
!!$
!!$             ponderation = ponderation + temp

          !Carré milieu
          !elseif ( ( u >= -1.0) .and. ( u <= 1.0) ) then
          if (u<=1) then

             temp = ((2-v)**3/6 - 4*(1-v)**3/6)*((2-u)**3/6 - 4*(1-u)**3/6)

             valeur = valeur + level(i)*temp

             ponderation = ponderation + temp
             compteur = compteur + 1
          !Carré droit
             !elseif ( ( u >= 1.0) .and. ( u <= 2.0) ) then
          elseif (abs(u-1.5)<=0.5) then

             temp = ((2-v)**3/6 - 4*(1-v)**3/6)*((2-u)**3/6)

             valeur = valeur + level(i)*temp

             ponderation = ponderation + temp
             compteur = compteur + 1

          end if


       !Ligne du bas
       !elseif ( ( v >= 1.0 ) .and. ( v <= 2.0 ) ) then
       elseif (abs(v-1.5)<=0.5) then

!!$          !Carré gauche
!!$          !if ( ( u >= -2.0) .and. ( u <= -1.0) ) then
!!$          if (abs(u+1.5)<=0.5) then
!!$
!!$             temp = ((2-v)**3/6)*((2-u)**3/6)
!!$
!!$             valeur = valeur + level(i)*temp
!!$
!!$             ponderation = ponderation + temp

          !Carré milieu
          !elseif ( ( u >= -1.0) .and. ( u <= 1.0) ) then
          if (u<=0.5) then

             temp = ((2-v)**3/6)*((2-u)**3/6 - 4*(1-u)**3/6)

             valeur = valeur + level(i)*temp

             ponderation = ponderation + temp
             compteur = compteur + 1

          !Carré droit
          !elseif ( ( u >= 1.0) .and. ( u <= 2.0) ) then
          elseif (abs(u-1.5)<=0.5) then

             temp = ((2-v)**3/6)*((2-u)**3/6)

             valeur = valeur + level(i)*temp

             ponderation = ponderation + temp
             compteur = compteur + 1

          end if


       end if

       if(compteur==25) exit

    end do

    if (ponderation > 1d-8) then  
       valeur = valeur/ponderation
    end if



  end subroutine noyau_interp





end module interpmod
