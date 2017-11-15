module kernel_mod
  implicit none
  real :: pi
  parameter (pi = 3.14159265)
  
contains
  
  function kfunc(dist,smdist)
    !cubic spline kernel function
    real :: kfunc
    real,intent(in) :: dist,smdist
    real :: q,sigw
    ! sigw = 2.0/(3.0*smdist**1)
    sigw=10.0/(7.0*pi*smdist**2)
    ! sigw = 1.0/(pi*smdist**3)
    q=dist/smdist
    if (q .lt. 1.0) then
       kfunc=1.0-(3.0/2.0)*(q)**2.0 + (3.0/4.0)*(q)**3.0
    else if (q .lt. 2.0) then
       kfunc = 0.25*(2.0-q)**3.0
    else
       kfunc = 0.0
    end if
    kfunc = sigw*kfunc
    return
  end function kfunc
  
  function dkfunc(dist,smdist)
    !cubic spline kernel function
    real :: dkfunc
    real,intent(in) :: dist,smdist
    real :: q,sigw
    !sigw=2.0/(3.0*smdist**1)
    sigw=10.0/(7.0*pi*smdist**2)
    !sigw=1.0/(pi*smdist**3)
    q = dist/smdist
    if (q .lt. 1.0) then
       dkfunc = -3.0*q + (9.0/4.0)*q**2.0
    else if (q .lt. 2.0) then
       dkfunc = -0.75*(2.0-q)**2.0
    else
       dkfunc = 0.0
    end if
    !	returns dw/dr. If dw/dh is needed,
    !multiply by (-r/h) outside
    dkfunc = sigw*dkfunc/smdist
    return
  end function dkfunc

!  function kfunc(dist,smdist)
!    ! Wendland kernel function
!    real :: kfunc
!    real,intent(in) :: dist,smdist
!    real :: q,sigw
!    sigw = 7.0/(4.0*pi*smdist**2.0)
!    q = dist/smdist
!    if(q.le.2.0)then
!       kfunc = sigw*(2.0*q + 1)*(1 - 0.5*q)**4.0
!    else
!       kfunc = 0.0
!    end if
!  end function kfunc
!
!  function dkfunc(dist,smdist)
!    ! Wendland kernel function derivative
!    real :: dkfunc
!    real, intent(in) :: dist,smdist
!    real q,sigw
!    sigw = 7.0/(4.0*pi*smdist**2.0)
!    if(q.le.2.0)then
!       dkfunc = (sigw*2.0/smdist)*((1 - 0.5*q)**4.0 - (2*q + 1)*(1 - 0.5*q)**3.0)
!    else
!       dkfunc = 0.0
!    end if
!  end function dkfunc

end module kernel_mod
