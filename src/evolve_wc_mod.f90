module evolve_wc_mod
  use particles_mod
  use neighbours_mod
  use bound_mod
  use kernel_mod
  use viscosity_mod
  implicit none
  
contains
  
  subroutine evolve_wcsph(nt)
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
    do i=1,n
       r0(i,:) = r(i,:)
       v0(i,:) = v(i,:)
       ro0(i) = ro(i)
    end do
    !set the time step
    call set_dt_wcsph
    
    ! predictor step
    do i=1,n
       v(i,:) = v0(i,:) + 0.5*dt*dvdt(i,:)
       r(i,:) = r0(i,:) + 0.5*dt*drdt(i,:)
       ro(i) = ro0(i) + 0.5*dt*drodt(i)
    end do
    
    ! find derivatives
    call calc_ddt_wcsph

    ! corrector step
    do i=1,n
       v(i,:) = v0(i,:) + 0.5*dt*dvdt(i,:)
       r(i,:) = r0(i,:) + 0.5*dt*drdt(i,:)
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
  end subroutine evolve_wcsph

  subroutine calc_ddt_wcsph
    real,dimension(dims) :: diss_v  !,f_st
    real,dimension(dims,dims) :: kg_cf
    real :: w,dw,vdotr,cij,roij,av
    real :: fij,fij0,rfij
    integer :: i,j,k,l

!    real :: st,s_ll,q_ij,F_int

    ! calculate tau if it's non-newtonian
    if(nnf) then
       call calculate_tau_nnf
    end if

    ! setup for tensile correction
    if(t_corr) then
       fij0 = 1.0/kfunc(h/h_ov_spacing,h)
    end if

    do i=1,n !for all internal particles
       drodt(i)=0.0
       dvdt(i,:) = 0.0 + g(:)
       drdt(i,:) = v(i,:)

       ! Calculate the kernel gradient correction terms, or don't
       if(kgc)then
          call kernel_gradient_correction(i,kg_cf)
       else
          ! identity matrix
          kg_cf(:,:) = 0.0
          do l=1,dims
             kg_cf(l,l) = 1.0
          end do
       end if
       
       do k=1,nneigh(i) !loop over neighbours
          j=neigh_list(i,k)

          ! calc distance, vdotr and 'ij' properties
          rij(:) = r(i,:) - r(j,:)
          vij(:) = v(i,:) - v(j,:)
          mrij = sqrt(dot_product(rij,rij))
          vdotr = dot_product(rij,vij)
          cij=0.5*(c(i)+c(j))
          roij=0.5*(ro(i)+ro(j))

          ! xsph adjustment
          if(xsph.ne.0.0)then
             w = kfunc(mrij,h)
             drdt(i,:) = drdt(i,:) + xsph*m(i)*vij(:)*w/roij
          end if

          ! tensile correction
          rfij = 0.0
          if(t_corr)then
             fij = kfunc(mrij,h)*fij0
             rfij = 0.2*(abs(min(p(i),0.0))/(ro(i)**2) + abs(min(p(j),0.0))/(ro(j)**2))
             rfij = rfij*fij**4.0
          end if

          ! find the gradient of kernel function
          dw = dkfunc(mrij,h)
          gradw(:) = dw*rij(:)/max(mrij,1e-6)
          gradw(:) = matmul(kg_cf,gradw)


          ! laminar viscosity or something...
          !diss_v(:) = 2.0*nu*m(i)*dot_product(rij,gradw)*vij(:)/(roij*(mrij**2 + 0.01*h**2))
          
          ! calculate the dissipative terms if it's a Non-Newtonian fluid
          if(nnf) then
             diss_v(:) = m(j)*matmul((tau(i,:,:)/ro(i)**2)+(tau(j,:,:)/ro(j)**2),gradw(:))
          else
             diss_v(:) = 0.0
          end if

          ! surface tension
          !st = 0.0728
          !s_ll = st/(0.0476*h_ov_spacing**4)
          !q_ij = mrij/h
          !F_int = s_ll*cos(3.0*pi*q_ij/4.0)
          !f_st(:) = F_int*rij(:)/(m(i)*(max(mrij,1e-6)))

          ! calculate AV terms if using AV.
          if(nnf.or.(vdotr.gt.0.0))then
             av = 0.0
          else
             av = av_value*cij*h*vdotr/(roij*(mrij**2.0+0.01*h**2.0))
          end if

          ! augment dvdt and drodt
          dvdt(i,:)=dvdt(i,:) - m(j)*gradw(:)* &
               (p(i)/(ro(i)**2.0) + p(j)/(ro(j)**2.0) + rfij + av)   + diss_v(:) ! + f_st(:)
          drodt(i) = drodt(i) + m(j)*vdotr*dw/max(mrij,1e-6)
       end do
    end do
    return
  end subroutine calc_ddt_wcsph
  
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
          w_cf = w_cf + w*m(j)/ro(j)  !kernel correction
          ro_new(i) = ro_new(i) + w*m(j)
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
    ! this can give some big errors in pressure...
    c0 = 20
    gm = 7.0
    B = c0*c0*ro_init/gm
           
    do i=1,n+ng
       p(i)= B*((ro(i)/ro_init)**gm - 1.0)
       c(i)=sqrt((c0**2)*(ro(i)/ro_init)**(gm-1.0))
    end do
    
    return
  end subroutine calc_p_from_ro

  subroutine set_dt_wcsph
    real vsig,minsf,minsc,tmp1
    real vdotr,acc
    integer i,j,k,l
    minsc=1e12
    minsf=1e12
    !loop over all internal particles
    do i=1,n
       do k=1,nneigh(i)
          j=neigh_list(i,k)

          ! courant constraint
          rij(:) = r(i,:) - r(j,:)
          vij(:) = v(i,:) - v(j,:)
          vdotr = dot_product(rij,vij)
          vsig=0.5*(c(i)+c(j)-1.0*vdotr)
          if(abs(h/vsig).le.minsc)then
             minsc=abs(h/vsig)
          end if
          !forces constraint
          acc=sqrt(dot_product(dvdt(i,:),dvdt(i,:)))
          if(abs(h/acc).le.minsf)then
             minsf=abs(h/acc)
          end if
       end do
    end do
    ! choose the most stringent
    dt=0.3*min(minsc,minsf)
    !write(6,*) "dtc",0.4*minsc,"dtf",0.25*minsf
    return
  end subroutine set_dt_wcsph

end module evolve_wc_mod
