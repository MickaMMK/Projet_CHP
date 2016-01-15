program main

  use transportmod
  use modmainmod
  use projection_methodmod
  use remplissage_poissonmod
  use raffinage

  implicit none

  !-----------------------------------!
  integer, parameter :: N = 50        !
  real(8), parameter :: tmax = 5      !
  real(8), parameter :: cfl = 0.9     !
  real(8), parameter :: period = 0.01 !
  !-----------------------------------!

  real(8), parameter :: pi = 3.1415926535897932384626433842795028841971693993

  real(8), dimension(N+1,N+1,2) :: noeuds, vitesses
  real(8), dimension(:,:), allocatable :: particules, vitesses_particules, new_part, partitemp
  real(8), dimension(N,N,2) :: centres
  real(8), dimension(N+1,N+1) :: level, rho, nu
  real(8), dimension(:), allocatable :: level_particules, lvlpartitemp
  real(8), dimension(N,N) :: pressions, rho_centre
  real(8), dimension(N*N,N*N) :: A
  integer, dimension(N*N) :: ipvt
  integer, dimension(:), allocatable :: npart_uni
  integer :: info
  real(8), dimension(2) :: g
  real(8) :: dx, dt, rho_air, rho_eau, nu_air, nu_eau, raff_size, lambda, t, tspent
  real(8) :: cfl_advection, cfl_visco, cfl_L2
  integer :: i, j, ci, di, ui, k, ki, meth, nbp, npart, raff_num, raff, temp, nbp_new, remaill, pos, transi, reproj
  character(len=1) :: c, d, u
  logical :: writo

  dx = 1./N
  npart = 0

  rho_air = 10. !800. !1.
  rho_eau = 1000. !1000.
  nu_air = 0.9d-2 !15.d-6
  nu_eau = 0.9d-2 !0.9d-6
  g = (/0.,-9.81/)

  print*, "Choix de la méthode :"
  print*, "1 - Eulerien"
  print*, "2 - Lagrangien"
  read*, meth
  nbp = 0
  if(meth < 1 .OR. meth > 2) then
     print*, "Erreur dans le choix de la méthode"
     print*, "==========="
     print*, "== ABORT =="
     print*, "==========="
     stop
  end if
  if(meth == 2) then
     print*, "Souhaitez-vous reprojeter les points lagrangiens sur la grille eulerienne"
     print*, "1 - Oui"
     print*, "2 - Non"
     read*, reproj
     if(reproj < 1 .OR. reproj > 2) then
        print*, "Erreur dans le choix de la reprojection"
        print*, "==========="
        print*, "== ABORT =="
        print*, "==========="
        stop
     else if(reproj == 2) then
        reproj = -1
     else if(reproj == 1) then
        print*, "Combien d'itérations entre chaque reprojection ?"
        read*, reproj
     end if
  end if
  print*, "Souhaitez-vous mettre l'eau :"
  print*, "1 - Au centre"
  print*, "2 - A l'extérieur"
  read*, pos
  if(pos < 1 .OR. pos > 2) then
     print*, "Erreur dans le choix de la position de l'eau"
     print*, "==========="
     print*, "== ABORT =="
     print*, "==========="
     stop
  elseif(pos == 1) then
     pos = -1
  elseif(pos == 2) then
     pos = 1
  end if
  print*, "Quelle transition à l'interface ?"
  print*, "1 - Discontinuité"
  print*, "2 - Linéaire"
  print*, "3 - Exponentielle"
  read*, transi
  if(transi < 1 .OR. transi > 3) then
     print*, "Erreur dans le choix de la transition à l'interface"
     print*, "==========="
     print*, "== ABORT =="
     print*, "==========="
     stop
  end if
  print*, "Souhaitez-vous rajouter des points lagrangiens aux abords de l'interface :"
  print*, "1 - Oui"
  print*, "2 - Non"
  read*, raff
  if(raff < 1 .OR. raff > 2) then
     print*, "Erreur dans le choix de la position de l'eau"
     print*, "==========="
     print*, "== ABORT =="
     print*, "==========="
     stop
  end if
     
  if(meth == 2) then
     nbp = (N-1)*(N-1)
  end if

  remaill = 2

  if(raff == 1) then
     print*, "Choisissez le nombre d'anneaux de chaque côté de l'interface :"
     read*, raff_num
     allocate(npart_uni(raff_num*2+1))
     print*, "Choisissez le nombre de points lagrangiens par anneau :"
     read*, npart_uni(1)
     npart_uni = npart_uni(1)
     print*, "Choisissez la distance entre chaque anneau :"
     read*, raff_size
     npart = npart_uni(1)*(raff_num*2+1)
     nbp = nbp + npart
     print*, "Voulez-vous remailler entre chaque itération ?"
     print*, "1 - Oui"
     print*, "2 - Non"
     read*, remaill
     if(remaill < 1 .OR. remaill > 2) then
        print*, "Erreur dans le choix du remaillage"
        print*, "==========="
        print*, "== ABORT =="
        print*, "==========="
        stop
     end if
  end if

  allocate(particules(nbp,2), vitesses_particules(nbp,2), level_particules(nbp))
  
  call initcoord(noeuds, centres)

  vitesses = 0.
  pressions = Patm
  
  do i = 1, N+1
     do j = 1, N+1
        level(i,j) = pos*(min(sqrt((noeuds(i,j,1)-0.5)**2+(noeuds(i,j,2)-0.5)**2) - 0.2, 10.d0))!noeuds(i,j,2)-0.1))
        if(sqrt((noeuds(i,j,1)-0.5)**2+(noeuds(i,j,2)-0.5)**2) - 0.2 .le. 0) then
           vitesses(i,j,2) = 0.
        end if
     end do
  end do

  if(meth == 2) then
     particules(1:(N-1)*(N-1),:) = vect2(noeuds(2:N,2:N,:))
     level_particules(1:(N-1)*(N-1)) = vect1(level(2:N,2:N))
  end if
  if(raff == 1) then
     do i = nbp-npart+1, nbp-npart+npart_uni(1)
        particules(i,:) = (/0.5+0.2*cos(2*pi*i/npart_uni(1)), 0.5+0.2*sin(2*pi*i/npart_uni(1))/)
        level_particules(i) = 0.
        do j = 1, raff_num
           particules(2*j*npart_uni(1)+i,:) = (/0.5+(0.2-j*raff_size)*cos(2*pi*i/npart_uni(1)), &
                & 0.5+(0.2-j*raff_size)*sin(2*pi*i/npart_uni(1))/)
           level_particules(2*j*npart_uni(1)+i) = -pos*raff_size*j
           particules((2*j-1)*npart_uni(1)+i,:) = (/0.5+(0.2+j*raff_size)*cos(2*pi*i/npart_uni(1)), &
                & 0.5+(0.2+j*raff_size)*sin(2*pi*i/npart_uni(1))/)
           level_particules((2*j-1)*npart_uni(1)+i) = pos*raff_size*j
        end do
     end do
  end if

  call write("level", 0, vect2(noeuds), vect1(level))
!!$  call write("vitesses_x", 0, vect2(noeuds), vect1(vitesses(:,:,1)))
!!$  call write("vitesses_y", 0, vect2(noeuds), vect1(vitesses(:,:,2)))
!!$  call write("pressions", 0, vect2(centres), vect1(pressions))

  if(raff == 1 .or. meth == 2) then
     call write("level_particules", 0, particules, level_particules)
!!$     call write("vitesses_x_particules", 0, particules, vitesses_particules(:,1))
!!$     call write("vitesses_y_particules", 0, particules, vitesses_particules(:,2))
  end if     

  t = 0.
  k = 0
  ki = 0
  tspent = 0
  writo = .false.

  lambda = 3 * dx

  cfl_visco = dx*dx/(4*max(nu_air,nu_eau))

!######################################################################################################################################################################################################
!######################################################################################################################################################################################################
!######################################################################################################################################################################################################

  do while (t .le. tmax)

     ki = ki + 1
!!$     if(ki == 3) then
!!$        stop
!!$     end if
     
     cfl_advection = dx/maxval(abs(vitesses))
     cfl_L2 = 2*min(nu_air,nu_eau)/maxval(abs(vitesses))**2
     dt = cfl*min(cfl_advection,cfl_visco,cfl_L2)

     if (tspent + dt .gt. period) then
        dt = period - tspent
        tspent = 0
        writo = .true.
        k = k + 1
     else
        if (period - (tspent + dt) .lt. 0.2d0*dt ) then
           dt = 0.7*dt
        end if
        tspent = tspent + dt
        writo = .false.
     end if

     print*, "==============================="
     print*, "dt = ",dt
     print*, "==============================="
  
     if(transi == 1) then

        rho = 0.5*sign(rho_air-rho_eau, level) + 0.5*(rho_air+rho_eau)
        nu = 0.5*sign(nu_air-nu_eau, level) + 0.5*(nu_air+nu_eau)

     ! transition linéaire

     else if(transi == 2) then

        do i = 1, N+1
           do j = 1, N+1

              rho(i,j) = min(rho_eau, max(rho_air, 0.5*(rho_air+rho_eau) + (1./lambda)*level(i,j)*0.5*(rho_air-rho_eau) )) 
              nu(i,j) = min(nu_eau, max(nu_air, 0.5*(nu_air+nu_eau) + (1./lambda)*level(i,j)*0.5*(nu_air-nu_eau) )) 

           end do
        end do

     ! transition "thick interface"

     else if(transi == 3) then

        do i = 1, N+1
           do j = 1, N+1

              rho(i,j) = min(max(rho_air,rho_eau), max(min(rho_air,rho_eau),  &
                   & min(rho_air,rho_eau)*(max(rho_air,rho_eau)/min(rho_air,rho_eau))**((level(i,j)+lambda)/(2*lambda))))                    
              nu(i,j) = min(max(nu_air,nu_eau), max(min(nu_air,nu_eau),  &
                   & min(nu_air,nu_eau)*(max(nu_air,nu_eau)/min(nu_air,nu_eau))**((level(i,j)+lambda)/(2*lambda))))                    

           end do
        end do

     end if


     rho_centre = (rho(1:N,1:N)+rho(1:N,2:N+1)+rho(2:N+1,1:N)+rho(2:N+1,2:N+1))/4

     if(meth == 1) then

        !==============================================================================================
        !========================================== EULERIEN ==========================================
        
        print*, 'Projection'
!!$        call projection_method(vitesses, pressions, rho, rho_centre, nu, g, dt, dx, level, A, ipvt)
        call projection_method_diphasique(vitesses, pressions, rho, rho_centre, rho*nu, g, dt, dx, level, A, ipvt)

        print*, 'Interpolation de la vitesse sur les particules'
        do i = 1, nbp
           call noyau_interp(vect2(noeuds), vect1(vitesses(:,:,1)), particules(i,:), dx, vitesses_particules(i,1))
           call noyau_interp(vect2(noeuds), vect1(vitesses(:,:,2)), particules(i,:), dx, vitesses_particules(i,2))
        end do

        print*, 'Transport'
        if(raff == 1) then
           call transport_level_EL(noeuds, vitesses, level, dt, dx, particules, level_particules)
           call transport_particules(particules, vitesses_particules, dt, dx)
           print*, 'Remaillage'
           call remaillage_particules(particules,vitesses_particules,level_particules,raff_num,nbp,npart,npart_uni,0.01d0,0.0001d0)
        else
           call transport_level(noeuds, vitesses, level, dt, dx)
        end if

        if(writo) then
           print*, 'Ecriture'
           call write("level", k, vect2(noeuds), vect1(level))
!!$           call write("vitesses_x", k, vect2(noeuds), vect1(vitesses(:,:,1)))
!!$           call write("vitesses_y", k, vect2(noeuds), vect1(vitesses(:,:,2)))
!!$           call write("pressions", k, vect2(centres), vect1(pressions))

           if(raff == 1) then
              call write("level_particules", k, particules, level_particules)
!!$              call write("vitesses_x_particules", k, particules, vitesses_particules(:,1))
!!$              call write("vitesses_y_particules", k, particules, vitesses_particules(:,2))
           end if
        end if

        !==============================================================================================

     else if(meth == 2) then
        
        !==============================================================================================
        !========================================= LAGRANGIEN =========================================
        
        print*, 'Méthode de projection'
!!$        call projection_method(vitesses, pressions, rho, rho_centre, nu, g, dt, dx, level, A, ipvt)
        call projection_method_diphasique(vitesses, pressions, rho, rho_centre, rho*nu, g, dt, dx, level, A, ipvt)

        print*, 'Interpolation de la vitesse sur les particules'
        do i = 1, nbp
           call noyau_interp(vect2(noeuds), vect1(vitesses(:,:,1)), particules(i,:), dx, vitesses_particules(i,1))
           call noyau_interp(vect2(noeuds), vect1(vitesses(:,:,2)), particules(i,:), dx, vitesses_particules(i,2))
        end do

        print*, 'Transport des particules'
        call transport_particules(particules, vitesses_particules, dt, dx)
        if(remaill == 1) then
           print*, 'Remaillage'
           call remaillage_particules(particules,vitesses_particules,level_particules,raff_num,nbp,npart,npart_uni,0.01d0,0.0001d0)
        end if

        print*, 'Interpolation du level sur les noeuds eulerien'
        do i = 2, N
           do j = 2, N
              call noyau_interp(particules, level_particules, noeuds(i,j,:), dx, level(i,j))
           end do
        end do

        if(reproj > 0) then
           if(mod(ki,reproj) == 0) then
              print*, 'Retour des points lagrangiens sur la grille eulerienne'
              particules(1:(N-1)*(N-1),:) = vect2(noeuds(2:N,2:N,:))
              vitesses_particules(1:(N-1)*(N-1),:) = vect2(vitesses(2:N,2:N,:))
              level_particules(1:(N-1)*(N-1)) = vect1(level(2:N,2:N))
           end if
        end if
        
        if(writo) then
           print*, 'Ecriture'
           call write("level_particules", k, particules, level_particules)
!!$           call write("vitesses_x_particules", k, particules, vitesses_particules(:,1))
!!$           call write("vitesses_y_particules", k, particules, vitesses_particules(:,2))

           call write("level", k, vect2(noeuds), vect1(level))
!!$           call write("vitesses_x", k, vect2(noeuds), vect1(vitesses(:,:,1)))
!!$           call write("vitesses_y", k, vect2(noeuds), vect1(vitesses(:,:,2)))
!!$           call write("pressions", k, vect2(centres), vect1(pressions))
        end if

        !==============================================================================================

     else
        
        print*, "Putain mais comment t'es arrivé là avec une mauvaise méthode ?!"

     end if

     t = t + dt

  end do

end program main
