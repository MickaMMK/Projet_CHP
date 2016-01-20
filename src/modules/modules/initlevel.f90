module initlevelmod

  implicit none

contains

!!$  subroutine initlevel
!!$
!!$    
!!$    
!!$    do i = 1, N+1
!!$       do j = 1, N+1
!!$          level(i,j) = pos*(sqrt((noeuds(i,j,1)-0.5)**2+(noeuds(i,j,2)-0.5)**2) - 0.2)
!!$          if(sqrt((noeuds(i,j,1)-0.5)**2+(noeuds(i,j,2)-0.5)**2) - 0.2 .le. 0) then
!!$             vitesses(i,j,2) = 0.
!!$          end if
!!$       end do
!!$    end do
!!$
!!$  end subroutine initlevel

end module initlevelmod
