module dtadapmod

  implicit none

contains
  
  subroutine dtadap(dt, tspent, period, writo, min, coeff, k)

    implicit none

    real(8), intent(inout) :: dt, tspent
    real(8), intent(in) :: period, min, coeff
    logical, intent(inout) :: writo    
    integer, intent(inout) :: k

    if (tspent + dt .gt. period) then
       dt = period - tspent
       tspent = 0
       writo = .true.
       k = k + 1
    else
       if (period - (tspent + dt) .lt. min*dt ) then
          dt = coeff*dt
       end if
       tspent = tspent + dt
       writo = .false.
    end if


  end subroutine dtadap

end module dtadapmod
