module grad_conj
  use matvectmod
  implicit none

  integer, parameter :: wp = kind(1.0D0)

contains

  subroutine grad_conj_opt(X,B,dx)
    implicit none

!    real(wp), dimension(:,:), intent(in)                  :: A
    real(wp), dimension(:), intent(in)                    :: B
    real(wp), dimension(size(B,1)), intent(out)           :: X
    real(wp), dimension(:,:), allocatable                 :: C, A_modif
    real(wp), dimension(:), allocatable                   :: w, r, Aw, B_modif
    real(wp), intent(in)                                  :: dx
    real(wp)                                              :: alpha, beta, eps, err
    integer                                               :: i, j, n

!    n = size(A,1)

   ! allocate(C(n,n), w(n), r(n), Aw(n))
    allocate(w(n), r(n), Aw(n))

!!$    C = 0._wp
!!$
!!$
!!$    do i = 1, n
!!$       C(i,i) = 1._wp / A(i,i) ! matmul(A,B) = mat_vect(B)
!!$    enddo
!!$
!!$    A_modif = matmul(C,A)
!!$    B_modif = matmul(C,B)

    eps = 0.000001_wp
    err = 1._wp

    X = 1._wp
    r = B - condi(mat_vect(X,dx))
    err = sqrt(dot_product(r,r))
    w = r
    alpha = dot_product(w,r) / dot_product(condi(mat_vect(w,dx)),r)
    i = 0

    do while((eps<err).and.(i<10000))

       i = i + 1
       X = X + alpha*w
       r = B - condi(mat_vect(X,dx))
       Aw = condi(mat_vect(w,dx))
       err = sqrt(dot_product(r,r))
       beta = dot_product(Aw-w,r) / dot_product(Aw,w)
       w = r - beta*w
       alpha = dot_product(w,r) / dot_product(condi(mat_vect(w,dx)),w)

    enddo

    print*,i

  end subroutine grad_conj_opt

  subroutine grad_conj_diphasique(X,B,dx,rho)
    implicit none

    real(wp), dimension(:), intent(in)                    :: B, rho
    real(wp), dimension(size(B,1)), intent(out)           :: X
    real(wp), dimension(:), allocatable                   :: h, g, g_old
!!$    real(wp), dimension(:), allocatable                   :: w, r, Aw
    real(wp), intent(in)                                  :: dx
    real(wp)                                              :: alpha, gamma, eps, err
    integer                                               :: i, j, n

    n = size(B)

    allocate(g(n), h(n), g_old(n))
!!$    allocate(w(n), r(n), Aw(n))

    eps = 0.0001_wp

    X = 1._wp

!!$    open(unit = 11, file = "matriceA", status = "replace")
!!$    open(unit = 12, file = "vectB", status = "replace")
!!$    do i = 1, size(B)
!!$       X = 0.d0
!!$       X(i) = 1.d0
!!$       write(11,*) mat_vect_diphasique(X,dx,rho)
!!$       write(12,*) B(i)
!!$    end do
!!$    close(11)
!!$    close(12)
!!$    stop


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Gradient conjugué avec conditionnement Pb: Pas symétrique. S'écrase mais diverge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!$    r = B - condi_diphasique(mat_vect_diphasique(X,dx,rho),rho)
!!$    err = sqrt(dot_product(r,r))
!!$    w = r
!!$    alpha = dot_product(w,r) / dot_product(condi_diphasique(mat_vect_diphasique(w,dx,rho),rho),r)
!!$    i = 0
!!$
!!$    do while((eps<err).and.(i<100000))
!!$
!!$       i = i + 1
!!$       X = X + alpha*w
!!$       r = B - condi_diphasique(mat_vect_diphasique(X,dx,rho),rho)
!!$       Aw = condi_diphasique(mat_vect_diphasique(w,dx,rho),rho)
!!$       err = sqrt(dot_product(r,r))
!!$       gamma = dot_product(Aw-w,r) / dot_product(Aw,w)
!!$       w = r - gamma*w
!!$       alpha = dot_product(w,r) / dot_product(condi_diphasique(mat_vect_diphasique(w,dx,rho),rho),w)
!!$
!!$    enddo
!!$
!!$    print*,i




!!$    r = B - mat_vect_diphasique(X,dx,rho)
!!$    err = abs(dot_product(condi_diphasique(r,rho),r))
!!$    w = condi_diphasique(r,rho)
!!$    alpha = dot_product(w,r) / dot_product(mat_vect_diphasique(w,dx,rho),r)
!!$    i = 0
!!$
!!$    do while((eps<err).and.(i<100000))
!!$
!!$       i = i + 1
!!$       X = X + alpha*w
!!$       r = B - mat_vect_diphasique(X,dx,rho)
!!$       Aw = condi_diphasique(mat_vect_diphasique(w,dx,rho),rho)
!!$       err = abs(dot_product(condi_diphasique(r,rho),r))
!!$       gamma = dot_product(Aw-w,r) / dot_product(Aw,w)
!!$       w = condi_diphasique(r,rho) - gamma*w
!!$       alpha = dot_product(w,r) / dot_product(mat_vect_diphasique(w,dx,rho),w)
!!$
!!$    enddo
!!$
!!$    print*,i
    
    g = mat_vect_diphasique(X,dx,rho) - B
    h = -condi_diphasique(g,rho)
    err = -dot_product(condi_diphasique(g,rho),g)
    i = 0

    do while((err > eps) .and. (i .le. 20000))

       alpha = -dot_product(g,h)/dot_product(h,mat_vect_diphasique(h,dx,rho))
       x = x + alpha*h
       g_old = g
       g = g + alpha*mat_vect_diphasique(h,dx,rho)
       gamma = dot_product(condi_diphasique(g,rho),g)/dot_product(condi_diphasique(g_old,rho),g_old)
       h = -condi_diphasique(g,rho) + gamma*h
       err = -dot_product(condi_diphasique(g,rho),g)
       i = i + 1

    end do

    print*, i


  end subroutine grad_conj_diphasique



!!$subroutine grad_conj_test(A,X,B)
!!$    implicit none
!!$
!!$   real(wp), dimension(:,:), intent(in)                  :: A
!!$    real(wp), dimension(:), intent(in)                    :: B
!!$    real(wp), dimension(size(A,1)), intent(out)           :: X
!!$    real(wp), dimension(:,:), allocatable                 :: Id
!!$    real(wp), dimension(:), allocatable                   :: g1, g2, d
!!$    real(wp)                                              :: alpha, beta, eps, err, w1
!!$    integer                                               :: i, j, n
!!$
!!$    n = size(A,1)
!!$
!!$    allocate(Id(n,n), g1(n), g2(n), d(n))
!!$
!!$    Id = 0._wp
!!$    do i = 1, n
!!$       Id(i,i) = 1._wp
!!$    enddo
!!$
!!$    w1 = 0.95_wp
!!$
!!$
!!$    X = -B
!!$    g1 = matmul(A,X) - B
!!$    d = g1
!!$
!!$    err = 1.0_wp
!!$
!!$    
!!$    do while(err>0.0001_wp)
!!$       alpha = 0.001_wp
!!$    !   do while(dot_product(d,(matmul((1._wp-w1)*A+Id),X)+alpha*d+(w1-1._wp)*B)>0._wp)
!!$    !      alpha = alpha / 2._wp
!!$    !   enddo
!!$
!!$       X = X + alpha*d
!!$       g2 = matmul(A,X) - B
!!$       beta = dot_product(g2,g2-g1) / dot_product(g1,g1)
!!$       d = -g2 + beta*d
!!$       g1 = g2
!!$
!!$       err = sqrt(dot_product(g2,g2))
!!$
!!$       print*,err
!!$
!!$    enddo
!!$
!!$
!!$
!!$  end subroutine grad_conj_test




end module grad_conj


