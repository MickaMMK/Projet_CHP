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

end module maillage