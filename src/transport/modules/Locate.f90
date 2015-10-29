module LOCATE
!Ce fichier est issu du TER de 1ere année de Mickaël et Pacôme. Il permet de trouver de manière rapide et efficace
!dans quel triangle se trouve un point du plan.
!On peut: -soit déclarer nos structures de données pour qu'elle soit identique à celle de ce fichier
!         -soit l'inverse.
  implicit none
  !%%%%%%%% DEFINITION DE TYPES %%%%%%%%!
  !---- Definition du maillage ----!
  ! type point
  type Vertex
     REAL, DIMENSION(2) :: X    ! Contient les coordonnees
  end type Vertex
     
  ! type triangle
  type Element
     INTEGER, DIMENSION(3) :: Sommet  ! Contient les sommets du triangle
     INTEGER, DIMENSION(3) :: Edge    ! Contient les aretes du triangle
     INTEGER, DIMENSION(3) :: Adj     ! Contient les triangles voisins
  end type Element

  ! type arete
  type Edge
     INTEGER, DIMENSION(2) :: Vertex ! Contient les noeuds
  end type Edge

  ! type maillage
  type Mesh2D
     INTEGER                                  :: Nv, Ne, Nedg  ! Nb de points, de triangles et d'aretes
     type(Vertex), DIMENSION(:), ALLOCATABLE  :: Vertex        ! Contient les sommets
     type(Element), DIMENSION(:), ALLOCATABLE :: Ele           ! Contient les triangles
     type(Edge), DIMENSION(:), ALLOCATABLE    :: Edge          ! Contient les aretes
  end type Mesh2D
  
  !---- Definition des hashtable pour la creation des voisinages ----!
  ! type item
  type Item2D 
     integer               :: Tri    ! 3*Num_tri + numero arete local
     integer, dimension(2) :: Vertex ! Numéro global des noeuds 
     integer               :: Nxt    ! pour dépacement
  end type Item2D
  
  ! type table de hash
  type HashTable2D
     type(Item2D), dimension(:), allocatable :: Item
     integer                                 :: SizeH ! Taille de la HashTable
     integer                                 :: SizeD ! Taille du dépacement
     integer                                 :: Compt ! Compteur de Nxt
  end type HashTable2D
  !%%%%%%%% FIN DE LA DECLARATION DES TYPES %%%%%%%%!
  
CONTAINS
  
  !==================================!
  SUBROUTINE ReadMesh(Mesh)
  !==================================!
 
    IMPLICIT NONE
    CHARACTER(len=99) :: fileName
    type(Mesh2D)      :: Mesh
    !------------------------------
    INTEGER  :: i, iter
    CHARACTER(len=99) :: MotCle
    type(HashTable2D) :: Hash
    !---------------------------

    OPEN ( unit = 13, file = 'parametre.data' )
    READ(13,*) FileName
    CLOSE(13)

    OPEN( unit = 12, file = TRIM(FileName)//'.mesh' ) 

    DO
       READ(12,*) MotCle
       IF ( TRIM(MotCle) == 'Vertices' ) THEN
          READ(12,*) Mesh%Nv
          READ(12,*) 
          ALLOCATE( Mesh%Vertex(Mesh%Nv) )
          DO i = 1, Mesh%Nv
             READ(12,*) Mesh%Vertex(i)%X(:)
          END DO
       ELSE IF ( TRIM(MotCle) == 'Triangles' ) THEN
          READ(12,*) Mesh%Ne
          READ(12,*)
          
          Hash%SizeH = 2*Mesh%Nv
          Hash%SizeD = 3*Mesh%Nv
          ALLOCATE ( Mesh%Ele(Mesh%Ne), Hash%Item(Hash%SizeH + Hash%SizeD) )
          
          DO iter = 1, Hash%SizeH + Hash%SizeD
             Hash%Item(iter)%Nxt = 0
             Hash%Item(iter)%Tri = -1
          END DO
          
          Hash%Compt = 1
          
          DO i = 1, Mesh%Ne
             READ(12,*) Mesh%Ele(i)%Sommet(:)
             Mesh%Ele(i)%Adj(1) = -1
             Mesh%Ele(i)%Adj(2) = -1
             Mesh%Ele(i)%Adj(3) = -1
             CALL FillHash2D(Mesh, i, Hash)
          END DO
          
       ELSE IF ( TRIM(MotCle) == 'Edges' ) THEN
          READ(12,*) Mesh%Nedg
          READ(12,*)
          ALLOCATE( Mesh%Edge(Mesh%Nedg) )
          DO i = 1, Mesh%Nedg
             READ(12,*) Mesh%Edge(i)%Vertex(:)
          END DO
       ELSE IF ( TRIM(MotCle) == 'End' ) THEN
          EXIT
          CLOSE(12)
       END IF
    END DO

  END SUBROUTINE ReadMesh
  !==================================!

  !=======================================!
  SUBROUTINE FillHash2D(Mesh, N_ele, Hash)
  !=======================================!
    
    IMPLICIT NONE
    type(Mesh2D)      :: Mesh
    type(HashTable2D) :: Hash
    INTEGER           :: N_ele
    !---------------------------

    INTEGER :: alpha, beta, cle, nd1, nd2, ModCle, pt1, pt2, NewLoc
    INTEGER :: i, j, iter, Next, Npts
    !-----------------------------------------------------------------

    alpha = 7
    beta = 11

    
    DO i = 1,2
       DO j = i+1,3
          pt1 = Mesh%Ele(N_ele)%Sommet(i)
          pt2 = Mesh%Ele(N_ele)%Sommet(j)
          IF (pt1 > pt2) THEN
             pt1 = Mesh%Ele(N_ele)%Sommet(j)
             pt2 = Mesh%Ele(N_ele)%Sommet(i)
          END IF

          cle = alpha*pt1 + beta*pt2
          ModCle = mod(cle,Hash%SizeH) + 1
          
          DO
             IF( Hash%Item(ModCle)%Tri == -1 ) THEN
                Hash%Item(ModCle)%Tri = 3*N_ele + 5 - (i + j)
                Hash%Item(ModCle)%Vertex(1) = pt1
                Hash%Item(ModCle)%Vertex(2) = pt2
                EXIT
             ELSE 
                nd1 = Hash%Item(ModCle)%Vertex(1)
                nd2 = Hash%Item(ModCle)%Vertex(2)
                IF ( (pt1 == nd1) .AND. (pt2 == nd2) ) THEN
                   Mesh%Ele(N_ele)%Adj(6 - (i + j)) = Hash%Item(ModCle)%Tri/3
                   Mesh%Ele(Hash%Item(ModCle)%Tri/3)%Adj(mod(Hash%Item(ModCle)%Tri,3) + 1) = N_ele
                   EXIT
                ELSE
                   IF ( Hash%Item(ModCle)%Nxt == 0 ) THEN
                      Hash%Item(ModCle)%Nxt = Hash%Compt
                      Hash%Compt = Hash%Compt + 1
                      ModCle = Hash%SizeH + Hash%Item(ModCle)%Nxt
                   ELSE
                      ModCle = Hash%SizeH + Hash%Item(ModCle)%Nxt
                   END IF
                END IF
             END IF
          END DO
       END DO
    END DO
    
  END SUBROUTINE FillHash2D
  !=======================================!

  !=========================================!
  SUBROUTINE LocateVertex2D(x, y, Mesh, Tri, BOOL)
  !=========================================!

    IMPLICIT NONE
    REAL         :: x, y
    type(Mesh2D) :: Mesh
    LOGICAL      :: BOOL
    
    INTEGER              :: Tri, iter
    REAL, DIMENSION(3)   :: CoorBar

    BOOL = .TRUE.

    Tri = 1
    iter = 1
    
    CALL CalcCoorBar(Tri, Mesh, x, y, CoorBar)

    DO WHILE ( ((CoorBar(1) < 0) .OR. (CoorBar(2) < 0) .OR. (CoorBar(3) < 0)) .AND. (iter <= Mesh%Ne) )
       
       IF (CoorBar(1) < 0) THEN
          IF (Mesh%Ele(Tri)%Adj(1) /= -1) THEN
             Tri = Mesh%Ele(Tri)%Adj(1)
          ELSE
             BOOL = .FALSE.
             EXIT
          END IF
       ELSE IF (CoorBar(2) < 0) THEN
          IF (Mesh%Ele(Tri)%Adj(2) /= -1) THEN
             Tri = Mesh%Ele(Tri)%Adj(2)
          ELSE
             BOOL = .FALSE.
             EXIT
          END IF
       ELSE
          IF (Mesh%Ele(Tri)%Adj(3) /= -1) THEN
             Tri = Mesh%Ele(Tri)%Adj(3)
          ELSE
             BOOL = .FALSE.
             EXIT
          END IF
       END IF
       
       CALL CalcCoorBar (Tri, Mesh, x, y, CoorBar)
       iter = iter + 1
    END DO
    
  END SUBROUTINE LocateVertex2D
  !=======================================!
  
  !=========================================!
  SUBROUTINE CalcCoorBar(Tri, Mesh, x, y, CoorBar)
  !=========================================!
    
    IMPLICIT NONE
    type(Mesh2D)         :: Mesh
    INTEGER              :: Tri
    REAL                 :: x, y
    REAL, DIMENSION(3)   :: CoorBar

    REAL    :: x1, y1, x2, y2, x3, y3, denom
    INTEGER :: S1, S2, S3
    
    S1 = Mesh%Ele(Tri)%Sommet(1)
    S2 = Mesh%Ele(Tri)%Sommet(2)
    S3 = Mesh%Ele(Tri)%Sommet(3)
    
    x1 = Mesh%Vertex(S1)%X(1)
    y1 = Mesh%Vertex(S1)%X(2)
    x2 = Mesh%Vertex(S2)%X(1)
    y2 = Mesh%Vertex(S2)%X(2)
    x3 = Mesh%Vertex(S3)%X(1)
    y3 = Mesh%Vertex(S3)%X(2)

    denom = (y - y3)*(x1 - x2)+(x - x3)*(y2 - y1)+(x - x2)*(y1 - y)+(y2 - y)*(x1 - x)
    
    CoorBar(1) = ( (y - y3)*(x - x2) + (x - x3)*(y2 - y) )/denom
    CoorBar(2) = ( (x3 - x)*(y1 - y) + (y - y3)*(x1 - x) )/denom
    CoorBar(3) = ( (x - x2)*(y1 - y) + (x1 - x)*(y2 - y) )/denom
    
  END SUBROUTINE CalcCoorBar
  !=========================================!





end program LOCATE
