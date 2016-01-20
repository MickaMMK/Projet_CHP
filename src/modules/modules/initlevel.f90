module initlevelmod

  implicit none

contains

  subroutine initlevel(N, noeuds, xc, yc, r, vcx, vcy, pos, level, vitesses)

    implicit none

    integer, intent(in) :: N, pos
    real(8), dimension(:,:,:), intent(in) :: noeuds
    real(8), intent(in) :: xc, yc, r, vcx, vcy
    real(8), dimension(:,:), intent(inout) :: level
    real(8), dimension(:,:,:), intent(inout) :: vitesses
    real(8) :: dist
    integer :: i, j
    
    do i = 1, N+1
       do j = 1, N+1
          dist = sqrt((noeuds(i,j,1)-xc)**2+(noeuds(i,j,2)-yc)**2) - r
          level(i,j) = pos*(dist)
          if(dist < 0) then
             vitesses(i,j,1) = vcx
             vitesses(i,j,2) = vcy
          end if
       end do
    end do

  end subroutine initlevel

end module initlevelmod
