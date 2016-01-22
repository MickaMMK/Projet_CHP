module chargemod
  implicit none

contains

   !---------------------------------------------------------------------------
   !> @brief
   !> Separate the tasks among the differents processor
   !
   !
   !> @param[in] me, n, Np    
   !> @param[out] i1, in   
   !> @return The number of the first and last task for each processor
   !--------------------------------------------------------------------------- 
 
  subroutine charge(me,n,Np,i1,in)

    implicit none
    !Variables d'entrées/sorties
    integer,intent(in) :: me, n, Np
    integer,intent(out) :: i1, in
    !Variables locales
    integer :: i,temp

    temp = n/Np

    i1 = 1
    in = i1 + temp-1
    do i=1,me
       temp = nint(float((n-in)/(Np-i)))
       i1 = in + 1
       in = i1 + temp-1
    end do

  end subroutine charge



end module chargemod
