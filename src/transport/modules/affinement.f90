module raffinage
  implicit none

contains
  subroutine remaillage(part,new_part,seuil_min,seuil_max)
    implicit none

    real(8), dimension(:,:), intent(in) :: part
    real(8), dimension(:,:),intent(out), allocatable :: new_part
    real(8), dimension(:,:), allocatable :: particules, new_particules
    real(8), dimension(:), allocatable :: distances, rajout
    real(8), intent(in) :: seuil_min, seuil_max
    real(8) :: dist, t1, t2
    integer :: i, j, k, N_stock, Npart_uni
    logical :: bool

    Npart_uni = size(part,1)
    allocate(new_particules(Npart_uni,2),particules(2*Npart_uni,2))
    particules(1:Npart_uni,:) = part
    particules(Npart_uni+1:2*Npart_uni,:) = part

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PREMIERE PARTIE, ON ENLEVE LES POINTS EN TROP !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    new_particules(1,:) = particules(1,:)
    i = 1
    k = 1
    bool = .FALSE.
    do while(i<Npart_uni)
       j = i + 1
       do while(sqrt((particules(i,1)-particules(j,1))**2+(particules(i,2)-particules(j,2))**2)<seuil_min)
          j = j + 1
          if(j>=Npart_uni) then
             bool = .TRUE.
             exit
          endif
       enddo
       if(bool) then
          exit
       endif
       i = j
       k = k + 1
       new_particules(k,:) = particules(i,:)
    enddo

    Npart_uni = k


    deallocate(particules)
    allocate(particules(Npart_uni,2),rajout(Npart_uni))
    particules(1:Npart_uni,:) = new_particules(1:Npart_uni,:)
    deallocate(new_particules)




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DEUXIEME PARTIE, ON REMPLIT LES TROUS !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    do while(k>0)

       rajout = 0.

       allocate(distances(Npart_uni))

       do i = 1, Npart_uni-1
          distances(i) = (particules(i+1,1)-particules(i,1))*(particules(i+1,1)-particules(i,1))+&
               &(particules(i+1,2)-particules(i,2))*(particules(i+1,2)-particules(i,2))
          if(distances(i)>seuil_max) then
             rajout(i) = 1
          endif
       enddo

       distances(Npart_uni) = (particules(1,1)-particules(Npart_uni,1))*(particules(1,1)-particules(Npart_uni,1))+&
            &(particules(1,2)-particules(Npart_uni,2))*(particules(1,2)-particules(Npart_uni,2))
       if(distances(Npart_uni)>seuil_max) then
          rajout(Npart_uni) = 1
       endif

       allocate(new_particules(Npart_uni+ceiling(sum(rajout)),2))

       j = 0

       do i = 1, Npart_uni-1
          new_particules(i+j,:) = particules(i,:)
          if(rajout(i)==1) then
             j = j + 1
             new_particules(i+j,:) = (particules(i,:) +  particules(i+1,:))*0.5
          endif
       enddo

       new_particules(Npart_uni+j,:) = particules(Npart_uni,:)
       if(rajout(i)==1) then
          j = j + 1
          new_particules(Npart_uni+j,:) = (particules(Npart_uni,:) +  particules(1,:))*0.5
       endif

       Npart_uni = Npart_uni+ceiling(sum(rajout))

       k = ceiling(sum(rajout))

       deallocate(particules,rajout,distances)
       allocate(particules(Npart_uni,2),rajout(Npart_uni))
       particules = new_particules
       deallocate(new_particules)

    enddo



!!!
    allocate(new_part(Npart_uni,2))

    new_part = particules



  end subroutine remaillage

  subroutine remaillage_particules(particules,vitesses_particules,level_particules,raff_num,nbp,npart,npart_uni,seuil1,seuil2)

    implicit none
    
    real(8), dimension(:,:), allocatable, intent(inout) :: particules, vitesses_particules
    real(8), dimension(:), allocatable, intent(inout) :: level_particules
    integer, dimension(:), intent(inout) :: npart_uni
    integer, intent(inout) :: raff_num, nbp, npart
    real(8), intent(in) :: seuil1, seuil2
    integer :: temp, nbp_new, i
    real(8), dimension(:,:), allocatable :: partitemp, new_part
    real(8), dimension(:), allocatable :: lvlpartitemp

    temp = 0
    do i = 1, raff_num*2+1
       call remaillage(particules(nbp-npart+temp+1:nbp-npart+temp+npart_uni(i),:), new_part, seuil1, seuil2)
       allocate(partitemp(nbp,2), lvlpartitemp(nbp))
       partitemp = particules
       lvlpartitemp = level_particules
       deallocate(particules, vitesses_particules, level_particules)
       allocate(particules(nbp-npart_uni(i)+size(new_part,1),2), vitesses_particules(nbp-npart_uni(i)+size(new_part,1),2))
       allocate(level_particules(nbp-npart_uni(i)+size(new_part,1)))
       nbp_new = nbp - npart_uni(i) + size(new_part,1)
       particules(1:nbp-npart+temp,:) = partitemp(1:nbp-npart+temp,:)
       level_particules(1:nbp-npart+temp) = lvlpartitemp(1:nbp-npart+temp)
       particules(nbp-npart+temp+1:nbp-npart+temp+size(new_part,1),:) = new_part
       level_particules(nbp-npart+temp+1:nbp-npart+temp+size(new_part,1)) = lvlpartitemp(nbp-npart+temp+1)
       particules(nbp-npart+temp+size(new_part,1)+1:nbp_new,:) = partitemp(nbp-npart+temp+npart_uni(i)+1:nbp,:)
       level_particules(nbp-npart+temp+size(new_part,1)+1:nbp_new) = lvlpartitemp(nbp-npart+temp+npart_uni(i)+1:nbp)
       nbp = nbp_new
       npart = npart - npart_uni(i) + size(new_part,1)
       npart_uni(i) = size(new_part,1)
       temp = temp + npart_uni(i)
       deallocate(new_part, partitemp, lvlpartitemp)
    end do

  end subroutine remaillage_particules

end module raffinage
