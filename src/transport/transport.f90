module transportmod
!!$%%%%%%% COMMENTAIRES %%%%%%%!!
!!$Résolution de l'équation de transport d(f)/dt + v.grad(f)=0
!!$
!!$ ATTENTION, CA DEVRAIT PAS COMPILER: VERIFIER SYNTAXE ET TYPES
!!$
!!$%%%%%%% FIN COMMENTAIRES %%%%%%%!!
  implicit none

contains
  subroutine transport(vitesses, level)
    real(8), dimension(N+1,N+1,2), intent(in) :: vitesses
    real(8), dimension(N+1,N+1), intent(inout) :: level
    integer::i,j,n
    type(maille),dimension(:),allocatable::f,f_old !contient la grandeur à transporter
    real(8)::t,dt

    !boucle sur les noeuds
    do i = 1, N+1
       do j = 1, N+1

          !position au temps précedent
          call euler(-1.*vitesses, point_old)
          call find_noeuds(point_old, indices)
          !valeur au temps précédent
          call interp(f_old%val(i),f%val(i), maille ou noeuds voisins)
          t=t+dt

       end do
    end do

  end subroutine transport

end module transportmod
