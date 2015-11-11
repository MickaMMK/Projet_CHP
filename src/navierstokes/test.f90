      program test
        use navierstokes
        implicit none

        real(kind=8) , dimension(9,9) :: A
        integer :: i


        call remplissage_poisson(A,1.d0,3)
        
        do i = 1 , 9

                write(*,"(9(F5.0,2X))"),A(i,:)


        end do

      end program test

