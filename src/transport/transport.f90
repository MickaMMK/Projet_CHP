module program transport
!!$%%%%%%% COMMENTAIRES %%%%%%%!!
!!$Résolution de l'équation de transport d(f)/dt + v.grad(f)=0
!!$
!!$ ATTENTION, CA DEVRAIT PAS COMPILER: VERIFIER SYNTAXE ET TYPES
!!$
!!$%%%%%%% FIN COMMENTAIRES %%%%%%%!!
  implicit none
integer::i,j,n
type(noeud),dimension(:),allocatable::vitesses
type(maille),dimension(:),allocatable::f,f_old !contient la grandeur à transporter
real(8)::t,dt

!boucle en temps
do while(t<tmax)

   f_old=f
   !boucle sur les noeuds
   do i=1,size(f_old%noeuds)
      !position au temps précedent
      call euler(f_old%noeud(i),vitesses)
      call locate(f_old%noeud, maille ou noeuds voisins)
      !valeur au temps précédent
      call interp(f_old%val(i),f%val(i), maille ou noeuds voisins)
t=t+dt

   end do
end do


end module transport
