!------------------------------------------------------------------------------
! ENSEIRB-MATMECA Colonne Splash 2015(c) - All rights reserved
!------------------------------------------------------------------------------

module remplissage_poissonmod
  use grad_conj
  implicit none

contains

   !---------------------------------------------------------------------------
   !> @brief
   !> Compute the \f$ N^2 \times N^2 \f$ Poisson's matrix with the pression condition on the last line.
   !
   !> @author 
   !> Corentin PRIGENT, engineer student at ENSEIRB-MATMECA, Bordeaux, FR.
   !
   !> @param[in] dx, N      
   !> @param[out] A      
   !> @return The Poisson's matrix
   !--------------------------------------------------------------------------- 
 
  !remplissage de la matrice du problème de Poisson
  subroutine remplissage_poisson(A,dx,N)
    implicit none

    real(kind=8) , dimension(:,:) , intent(out) :: A
    real(kind=8) , intent(in) :: dx
    integer , intent(in) :: N !matrice de pression de taille NxN

    integer :: i , j 

    A = 0

    !premier bloc de N lignes

    A(1,1) = -2
    A(1,2) = 1
    A(1,N+1) = 1

    do i = 2 , N-1 

       A(i,i) = -3
       A(i,i+1) = 1
       A(i,i-1) = 1
       A(i,i+N) = 1

    end do

    A(N,N) = -2
    A(N,N-1) = 1
    A(N,2*N) = 1

    !dernier bloc de N lignes

    A(N*(N-1)+1,N*(N-1)+1) = -2
    A(N*(N-1)+1,N*(N-1)+2) = 1
    A(N*(N-1)+1,N*(N-1)+1-N) = 1

    do i = 2 , N-1 

       A(N*(N-1)+i,N*(N-1)+i) = -3
       A(N*(N-1)+i,N*(N-1)+i+1) = 1
       A(N*(N-1)+i,N*(N-1)+i-1) = 1
       A(N*(N-1)+i,N*(N-1)+i-N) = 1

    end do

    A(N*N,N*N) = dx*dx

    !blocs génériques

    do i = 2 , N-1

       A((i-1)*N+1,(i-1)*N+1) = -3
       A((i-1)*N+1,(i-1)*N+2) = 1
       A((i-1)*N+1,(i-1)*N+1+N) = 1
       A((i-1)*N+1,(i-1)*N+1-N) = 1

       do j = 2 , N-1

          A((i-1)*N+j,(i-1)*N+j) = -4
          A((i-1)*N+j,(i-1)*N+j+1) = 1
          A((i-1)*N+j,(i-1)*N+j-1) = 1
          A((i-1)*N+j,(i-1)*N+j+N) = 1
          A((i-1)*N+j,(i-1)*N+j-N) = 1

       end do
       
       A(i*N,i*N) = -3
       A(i*N,i*N-1) = 1
       A(i*N,i*N+N) = 1
       A(i*N,i*N-N) = 1

    end do

  end subroutine remplissage_poisson

end module remplissage_poissonmod
