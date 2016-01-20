module initpartmod

  implicit none

  real(8), parameter :: pi = 3.1415926535897932384626433842795028841971693993

contains

  subroutine initpart(part_raff, level_part_raff, npart_uni, raff_num, raff_size, pos)

    real(8), dimension(:,:), intent(inout) :: part_raff
    real(8), dimension(:), intent(inout) :: level_part_raff
    integer, intent(in) :: npart_uni, raff_num, pos
    real(8), intent(in) :: raff_size
    integer :: i, j
    
    do i = 1, npart_uni
       part_raff(i,:) = (/0.5+0.2*cos(2*pi*i/npart_uni), 0.5+0.2*sin(2*pi*i/npart_uni)/)
       level_part_raff(i) = 0.
       do j = 1, raff_num
          part_raff(2*j*npart_uni+i,:) = (/0.5+(0.2-j*raff_size)*cos(2*pi*i/npart_uni), &
               & 0.5+(0.2-j*raff_size)*sin(2*pi*i/npart_uni)/)
          level_part_raff(2*j*npart_uni+i) = -pos*raff_size*j
          part_raff((2*j-1)*npart_uni+i,:) = (/0.5+(0.2+j*raff_size)*cos(2*pi*i/npart_uni), &
               & 0.5+(0.2+j*raff_size)*sin(2*pi*i/npart_uni)/)
          level_part_raff((2*j-1)*npart_uni+i) = pos*raff_size*j
       end do
    end do

  end subroutine initpart

end module initpartmod
