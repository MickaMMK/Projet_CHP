module chargemod
  implicit none

contains

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
