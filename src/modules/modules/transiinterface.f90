module transiinterfacemod

  implicit none

contains

  subroutine transiinterface(N, transi, level, rho, nu, rho_eau, rho_air, nu_eau, nu_air, lambda)
    
    implicit none

    integer, intent(in) :: N, transi
    real(8), dimension(:,:), intent(in) :: level
    real(8), intent(in) :: rho_eau, rho_air, nu_eau, nu_air, lambda
    real(8), dimension(:,:), intent(inout) :: rho, nu
    integer :: i, j

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

  end subroutine transiinterface

end module transiinterfacemod
