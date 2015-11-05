module maillage
  
  implicit none

  type noeud
     real(8), dimension(2) :: coord
     type(maille), dimension(:), pointer :: mailles
  end type noeud

  type centre
     real(8), dimension(2) :: coord
     type(maille), pointer :: maille
  end type centre

  type maille
     integer :: num
     type(noeud), dimension(:), pointer :: sommets
     type(centre), pointer :: centre
  end type maille

  type vector2D
        real(kind=8) , dimension(2) :: coord
        type(noeud) , pointer :: noeud
  end type vector2D

  

end module maillage
