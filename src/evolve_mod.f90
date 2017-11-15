module evolve_mod
  use particles_mod
  use neighbours_mod
  use bound_mod
  use kernel_mod
  implicit none
  
contains
  
  subroutine evolve(nt)
    ! this is the "leapfrog" algorithm, see Price PhD
    ! Predictor step
    ! v_1 = v_0 + 0.5*dt*dvdt_0
    ! Corrector step
    ! v_2 = v_0 + 0.5*dt*dvdt_1
    ! Final step
    ! v = 2*v_2 - v_0
    ! In practice, replace dvdt_0 with dvdt_1 from previous
    ! time-step, and ddt, neighbours, etc only
    ! need to be called once per time-step!
    integer,intent(in) :: nt
    real,allocatable,dimension(:) :: ro0
    real,allocatable,dimension(:,:) :: r0,v0
    integer :: i

    allocate(ro0(n))
    allocate(v0(n,dims))
    allocate(r0(n,dims))
    
    ! put variables in the "0" arrays
    r0(:,:) = r(1:n,:)
    v0(:,:) = v(1:n,:)
    ro0(:) = ro(1:n)
    
    ! predictor step
    do i=1,n
       v(i,:) = v0(i,:) + 0.5*dt*dvdt(i,:)
       r(i,:) = r0(i,:) + 0.5*dt*v(i,:)
       ro(i) = ro0(i) + 0.5*dt*drodt(i)
    end do
    
    ! find derivatives
    call calc_ddt

    ! corrector step
    do i=1,n
       v(i,:) = v0(i,:) + 0.5*dt*dvdt(i,:)
       r(i,:) = r0(i,:) + 0.5*dt*v(i,:)
       ro(i) = ro0(i) + 0.5*dt*drodt(i)
    end do

    !final step
    do i=1,n
       v(i,:) = 2.0*v(i,:) - v0(i,:)
       r(i,:) = 2.0*r(i,:) - r0(i,:)
       ro(i) = 2.0*ro(i) - ro0(i)
    end do
    
    deallocate(ro0)
    deallocate(v0)
    deallocate(r0)

    ! apply boundary conditions
    call boundaries
    
    ! find neighbours
    call find_neighbours
    
    ! density re-initialisation every 30 time steps
    if(mod(nt,30).eq.0)then
       call density_reinitialisation
    end if

    ! calculate the pressure
    call calc_p_from_ro
    
    return
  end subroutine evolve

  subroutine calc_ddt
    ! subroutine to calculate d/dt for v,ro, and also dt.
    real :: dw,vdotr,cij,roij,av
    real :: minsf,minsc,vsig,acc
    integer :: i,j,k

    ! setup for calculating dt
    minsc=1e12
    minsf=1e12

    do i=1,n !for all internal particles
       drodt(i) = 0.0
       dvdt(i,:) = 0.0 + g(:)
       
       do k=1,nneigh(i) !loop over neighbours
          j=neigh_list(i,k)

          ! calc distance, vdotr and 'ij' properties
          rij(:) = r(i,:) - r(j,:)
          vij(:) = v(i,:) - v(j,:)
          mrij = sqrt(dot_product(rij,rij))
          vdotr = dot_product(rij,vij)
          cij=0.5*(c(i)+c(j))
          roij=0.5*(ro(i)+ro(j))

          ! time-step calculation - Courant constraint
          vsig = cij - 0.5*vdotr
          if(abs(h/vsig).le.minsc)then
             minsc=abs(h/vsig)
          end if

          ! calculate AV terms if using AV.
          if(vdotr.gt.0.0)then
             av = 0.0
          else
             av = av_value*cij*h*vdotr/(roij*(mrij**2.0+0.01*h**2.0))
          end if
          
          ! find the gradient of kernel function
          dw = dkfunc(mrij,h)
          gradw(:) = dw*rij(:)/max(mrij,1e-6)

          ! augment dvdt and drodt
          dvdt(i,:) = dvdt(i,:) - m*gradw(:)*(p(i)/(ro(i)**2.0) + p(j)/(ro(j)**2.0)  + av)
          drodt(i) = drodt(i) + m*vdotr*dw/max(mrij,1e-6)
       end do
       
       ! time-step calculation - forces constraint
       acc=sqrt(dot_product(dvdt(i,:),dvdt(i,:)))
       if(abs(h/acc).le.minsf)then
          minsf=abs(h/acc)
       end if
       
    end do
    
    ! choose the most stringent time-step constraint
    dt=0.3*min(minsc,minsf)
    
    return
  end subroutine calc_ddt
  
  subroutine density_reinitialisation
    ! As in Gomez-Gesteira et al, 2012
    ! NB, uses a correction factor to the kernel to account for free surfaces
    real, allocatable,dimension(:) :: ro_new
    real :: w,w_cf
    integer :: i,j,k

    allocate(ro_new(n))
    
    do i=1,n
       w_cf = 0.0
       ro_new(i) = 0.0
       do k=1,nneigh(i) !loop over neighbours
          j=neigh_list(i,k)
          mrij = sqrt(dot_product(r(i,:)-r(j,:),r(i,:)-r(j,:)))
          w = kfunc(mrij,h)
          w_cf = w_cf + w*m/ro(j)  !kernel correction
          ro_new(i) = ro_new(i) + w*m
       end do
       ro_new(i) = ro_new(i)/w_cf
    end do
    ! pass ro_new back into first n elements of ro
    ro(:) = ro_new(:)

    deallocate(ro_new)
    return
  end subroutine density_reinitialisation

  subroutine calc_p_from_ro
    integer :: i
    real :: c0,B,gm

    ! calculate p and c, using Tait EoS
    c0 = 20 ! lower c0 gives higher dt and bigger errors in pressure
    gm = 7.0
    B = c0*c0*ro_init/gm
           
    do i=1,n
       p(i)= B*((ro(i)/ro_init)**gm - 1.0)
       c(i)=sqrt((c0**2)*(ro(i)/ro_init)**(gm-1.0))
    end do
    
    return
  end subroutine calc_p_from_ro

end module evolve_mod
