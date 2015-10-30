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

!  type(maille), target :: case
!  type(centre), target :: point
!  type(noeud), dimension(4), target :: sommet

!  case%num = 1
!  point%coord = (/0.5,0.5/)

!  case%centre => point
!  point%maille => case

!  print*, "Le centre de la case",case%num,"est",case%centre%coord

!  sommet(1)%coord = (/0,0/)
!  sommet(2)%coord = (/1,0/)
!  sommet(3)%coord = (/1,1/)
!  sommet(4)%coord = (/0,1/)

!  case%sommets => sommet

  

!  print*, "On prend le centre de coordonnées",point%coord,", il correspond au centre&
!       & de la case",point%maille%num," dont les sommets sont les suivants :"
!  print*, point%maille%sommets(1)%coord
!  print*, point%maille%sommets(2)%coord
!  print*, point%maille%sommets(3)%coord
!  print*, point%maille%sommets(4)%coord


end module maillage
