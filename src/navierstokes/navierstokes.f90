module navierstokes
     use maillage
     use vecteur_plan
     implicit none

     contains

!méthode de projection de Chorin
     subroutine projection_method(u,p,u_next,p_next)
        implicit none

        type(vector2D) , dimension(:,:) , intent(in) :: u
        real(kind=8) , dimension(:,:) , intent(in) :: p

        type(vector2D) , dimension(:,:), intent(out) :: u_next
        real(kind=8) , dimension(:,:) , intent(out) :: p_next

        type(vector2D) , dimension(:,:) :: u_star

        u_star = u 

     end subroutine projection_method

!schéma MAC
     subroutine MAC_scheme()
        implicit none


     end subroutine MAC_scheme

!schéma Adams-Bashforth + Crank-Nicolson
     subroutine troisieme_schema()
        implicit none


     end subroutine troisieme_schema

module navierstokes
