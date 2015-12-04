!------------------------------------------------------------------------------
! NASA/GSFC, Software Integration & Visualization Office, Code 610.3
!------------------------------------------------------------------------------
!
! MODULE: Module Name
!
!> @author
!> Module Author Name and Affiliation
!
! DESCRIPTION: 
!> Brief description of module.
!
! REVISION HISTORY:
! DD Mmm YYYY - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------

program main

  use transportmod
  use modmainmod
  use projection_methodmod
  use remplissage_poissonmod
  
  implicit none

  !--------------------------------
  integer, parameter :: N = 50     !
  integer, parameter :: tmax = 5  !
  real(8), parameter :: cfl = 0.9 !
  !--------------------------------

  real(8), dimension(N+1,N+1,2) :: noeuds, vitesses
  real(8), dimension(N,N,2) :: centres
  real(8), dimension(N+1,N+1) :: level, rho, nu
  real(8), dimension(N,N) :: pressions, rho_centre
  real(8), dimension(N*N,N*N) :: A
  integer, dimension(N*N) :: ipvt
  integer :: info
  real(8), dimension(2) :: g
  real(8) :: dx, dt, rho_air, rho_eau, nu_air, nu_eau
  integer :: i, j, Niter, ci, di, ui, k
  character(len=1) :: c, d, u
  
  dx = 1./N

  rho_air = 1000. !1.
  rho_eau = 1000. !1000.
  nu_air = 0.9E-2 !15-6
  nu_eau = 0.9E-2 !0.9E-6
  g = (/0.,-9.81/)
  
  call initcoord(noeuds, centres)

  vitesses = 0.
  vitesses(:,N+1,1) = 2.
!!$  vitesses(2:N,2:N,:) = 0.3

!!$  do i = 2, N
!!$     do j = 2, N
!!$        vitesses(i,j,:) = (/-(noeuds(i,j,1)-0.5)**2+0.25,-(noeuds(i,j,2)-0.5)**2+0.25/)
!!$        print*, vitesses(i,j,:)
!!$     end do
!!$  end do
  pressions = 1.013D5
  
  do i = 1, N+1
     do j = 1, N+1
        if(sqrt((noeuds(i,j,1)-0.25)**2+(noeuds(i,j,2)-0.75)**2) < 0.2) then
           level(i,j) = 0.
        else
           level(i,j) = 1.
        end if
     end do
  end do

  print*, maxval(abs(vitesses))
  !dt = cfl*dx/5
  dt = cfl*min(dx/maxval(abs(vitesses)),dx*dx/(4*max(nu_air,nu_eau)),2*min(nu_air,nu_eau)/maxval(abs(vitesses))**2)
  print*, "dt = ",dt
  Niter = ceiling(tmax/dt)

  call write("level", 0, noeuds, level)
  call write("vitesses_x", 0, noeuds, vitesses(:,:,1))
  call write("vitesses_y", 0, noeuds, vitesses(:,:,2))
  call write("pressions", 0, centres, pressions)

  call remplissage_poisson(A,dx,N)
  call DGETRF(N*N, N*N, A, N*N, ipvt, info)

  do k = 1, Niter

     print*, "Itération ",k," sur ",Niter
  
     rho = int(level+0.5)*rho_air + (1-int(level+0.5))*rho_eau
     nu = int(level+0.5)*nu_air + (1-int(level+0.5))*nu_eau
     rho_centre = (rho(1:N,1:N)+rho(1:N,2:N+1)+rho(2:N+1,1:N)+rho(2:N+1,2:N+1))/4
     !print*, pressions

     print*, 'Projection'
     call projection_method(vitesses, pressions, rho, rho_centre, nu, g, dt, dx, level,A,ipvt)

     print*, 'Transport'
     call transport(noeuds, vitesses, level, dt, dx)

     print*, 'Ecriture'
     call write("level", k, noeuds, level)
     call write("vitesses_x", k, noeuds, vitesses(:,:,1))
     call write("vitesses_y", k, noeuds, vitesses(:,:,2))
     call write("pressions", k, centres, pressions)

     !dt = cfl*min(dx/maxval(vitesses),dx*dx/(2*max(nu_air,nu_eau)))

  end do

contains

   !---------------------------------------------------------------------------  
   !> @author 
   !> Routine Author Name and Affiliation.
   !
   ! DESCRIPTION: 
   !> Brief description of routine. 
   !> @brief
   !> Flow method (rate of change of position) used by integrator.
   !> Compute \f$ \frac{d\lambda}{dt} , \frac{d\phi}{dt},  \frac{dz}{dt} \f$
   !
   ! REVISION HISTORY:
   ! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
   !
   !> @param[in] inParam      
   !> @param[out] outParam      
   !> @return returnValue
   !--------------------------------------------------------------------------- 

end program main
