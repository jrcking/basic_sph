!! this file is part of Jack King's basic SPH code for simulating a Dam Break
!! 

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
    sigw=10.0/(7.0*pi*smdist**2) ! normalisation for 2D hard-coded here!
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
    sigw=10.0/(7.0*pi*smdist**2) ! normalisation for 2D hard-coded here!
    q = dist/smdist
    if (q .lt. 1.0) then
       dkfunc = -3.0*q + (9.0/4.0)*q**2.0
    else if (q .lt. 2.0) then
       dkfunc = -0.75*(2.0-q)**2.0
    else
       dkfunc = 0.0
    end if
    ! This returns dw/dr. If dw/dh is needed (for example, when h varies
    ! with time), multiply by (-r/h) outside this subroutine.
    dkfunc = sigw*dkfunc/smdist
    return
  end function dkfunc
  
end module kernel_mod
