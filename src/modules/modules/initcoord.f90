module initcoordmod

  implicit none

contains

  subroutine initcoord(noeuds, centres)

    implicit none

    real(8), dimension(:,:,:), intent(inout) :: noeuds, centres
    real(8) :: dx
    integer :: i, j

    dx = 1./size(centres,1) 
    
    do i = 1, size(noeuds,1)
       do j = 1, size(noeuds,2)
          noeuds(i,j,:) = (/(i-1)*dx, (j-1)*dx/)
       end do
    end do

    do i = 1, size(centres,1)
       do j = 1, size(centres,2)
          centres(i,j,:) = (/((i-1)+0.5)*dx, ((j-1)+0.5)*dx/)
       end do
    end do

  end subroutine initcoord

end module initcoordmod
