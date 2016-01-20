module ordonnancement
  implicit none

  type ordpoint
     real(8) :: x, y
     integer :: numero
  end type ordpoint

contains

  subroutine reordre(particules)
    implicit none
    type(ordpoint), dimension(:), intent(inout) :: particules
    type(ordpoint), dimension(:), allocatable :: temp
    integer :: i, n

    n = size(particules)
    allocate(temp(n))

    do i = 1, n
       temp(particules(i)%numero)%x = particules(i)%x
       temp(particules(i)%numero)%y = particules(i)%y
       temp(i)%numero = i
    end do

    particules = temp
  end subroutine reordre

end module ordonnancement
