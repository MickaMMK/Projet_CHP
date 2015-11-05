module eulermod
!!$%%%%%%% COMMENTAIRES %%%%%%%!!
!!$ Schéma d'Euler pour situer au temps n-1 la valeur sur un noeud au temps n
!!$ de la fonction cherchée 
!!$
!!$ ATTENTION, CA DEVRAIT PAS COMPILER: VERIFIER SYNTAXE ET TYPES
!!$
!!$%%%%%%% FIN COMMENTAIRES %%%%%%%!!
  implicit none

contains
  subroutine euler(vitesses,f,dt)
    implicit none
    type(vecteur),dimension(:),allocatable::vitesses
    type(scalaire),dimension(:),allocatable::f,f_old !contient la grandeur à transporter
    integer::i,j
    real(8)::dt
    !on a besoin de la vitesse au temps n
    f_old=f

    f_old%noeuds(i)=f_old%noeuds(i)-vitesses(i)*dt
    !ca devrait donner la position au temps n, qui est le temps preced
  end subroutine euler


end module eulermod
