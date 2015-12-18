module vectmod

  implicit none

contains

  function vect1(array) result(vector)

    implicit none
    
    real(8), dimension(:,:), intent(in) :: array
    real(8), dimension(size(array,1)*size(array,2)) :: vector
    integer :: i, j

    do i = 1, size(array,1)
       do j = 1, size(array,2)
          vector(i+size(array,1)*(j-1)) = array(i,j)
       end do
    end do

  end function vect1

  function vect2(array) result(vector)

    implicit none
    
    real(8), dimension(:,:,:), intent(in) :: array
    real(8), dimension(size(array,1)*size(array,2),size(array,3)) :: vector
    integer :: i, j

    do i = 1, size(array,1)
       do j = 1, size(array,2)
          vector(i+size(array,1)*(j-1),:) = array(i,j,:)
       end do
    end do

  end function vect2

  function unvect1(vector) result(array)

    implicit none
    
    real(8), dimension(:), intent(in) :: vector
    real(8), dimension(int(sqrt(float(size(vector)))),int(sqrt(float(size(vector))))) :: array
    integer :: i, j

    do i = 1, size(array,1)
       do j = 1, size(array,2)
          array(i,j) = vector(i+size(array,1)*(j-1))
       end do
    end do

  end function unvect1

  function unvect2(vector) result(array)

    implicit none
    
    real(8), dimension(:,:), intent(in) :: vector
    real(8), dimension(int(sqrt(float(size(vector,1)))),int(sqrt(float(size(vector,1)))),size(vector,2)) :: array
    integer :: i, j

    do i = 1, size(array,1)
       do j = 1, size(array,2)
          array(i,j,:) = vector(i+size(array,1)*(j-1),:)
       end do
    end do

  end function unvect2

end module vectmod
