program navierstokes
        use maillage
        implicit none

        type(noeud) :: test
        test%coord = (/1.,2./)

        write(*,*)test%coord


end program navierstokes
