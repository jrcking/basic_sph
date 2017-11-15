module evolve_i_mod
  use particles_mod
  use neighbours_mod
  use bound_mod
  use kernel_mod
  use viscosity_mod
  implicit none
  

contains
  subroutine evolve_isph(nt)
    ! this is isph as described in Lind et al 2015
    integer,intent(in)::nt
    real,allocatable,dimension(:,:)::r0,v0
    integer i,l

    allocate(v0(n,dims))
    allocate(r0(n,dims))
    
    ! put variables in the "0" arrays
    do i=1,n
       r0(i,:) = r(i,:)
       v0(i,:) = v(i,:)
    end do

    !set the time step
    call set_dt_isph

    ! find the time-derivatives
    call calc_ddt_isph
    
    ! projection step
    do i=1,n
       r(i,:) = r0(i,:) + dt*v(i,:)
       v(i,:) = v0(i,:) + dt*dvdt(i,:)
    end do

    ! apply boundary conditions and find neighbours
    call boundaries
    call find_neighbours

    ! solve a PPE to obtain the pressure, and return dvdt=(1/ro)*grad(p)
    call ppe
    
    ! corrector step
    do i=1,n
       v(i,:) = v(i,:) + dt*dvdt(i,:)
       r(i,:) = r0(i,:) + 0.5*dt*(v(i,:) + v0(i,:))
    end do

    ! apply boundary conditions and find neighbours
    call boundaries
    call find_neighbours

    deallocate(v0)
    deallocate(r0) 

    ! fickian shifting...
    call fickian_shift
    
    ! apply boundary conditions and find neighbours
    call boundaries
    call find_neighbours
    
    return
  end subroutine evolve_isph

  subroutine ppe
    ! subroutine based on Xu PhD thesis amongst others.
    ! sets up the problem AP=B and solves for P.
    real,dimension(dims,dims):: kg_cf
    real,allocatable,dimension(:,:)::A
    real,allocatable,dimension(:)::B,X
    logical,dimension(n+ng)::fsp
    real dw,div_r,tmp
    integer i,j,k,l

    allocate(A(n+ng,n+ng))
    allocate(B(n+ng))
    allocate(X(n+ng))

    A(:,:) = 0.0
    B(:) = 0.0
 
    do i=1,n+ng
       call kernel_gradient_correction(i,kg_cf)
       div_r = 0.0
       do k=1,nneigh(i)
          j=neigh_list(i,k)

          ! find some ij-properties
          rij(:) = r(i,:) - r(j,:)
          vij(:) = v(i,:) - v(j,:)
          mrij = sqrt(dot_product(rij,rij))
          dw = dkfunc(mrij,h)
          gradw(:) = dw*rij(:)/max(mrij,1e-6)

          ! calculate div.r (is it free surface particle?)
          div_r = div_r - m(j)*dot_product(rij,gradw)/ro(j)

          ! correct kernel gradient for populating A,B
          gradw(:) = matmul(kg_cf,gradw)

          ! populate A
          A(i,j) = -2.0*m(j)*dot_product(rij,gradw)/ &
               (ro(j)*ro(j)*(mrij**2 + 0.01*h**2))

          ! populate B
          B(i) = B(i) + m(j)*dot_product(-1.0*vij,gradw)/(ro(j)*dt)
       end do
       ! is it FSP?
       if(div_r.le.1.5)then
          fsp(i) = .true.
       else
          fsp(i) = .false.
       end if
       A(i,i) = -1.0*sum(A(i,:))
    end do

    ! apply the free surface boundary condition..
    do i=1,n
       if(fsp(i))then
          A(i,:) = 0.0
          A(:,i) = 0.0
          B(i) = 0.0
       end if       
    end do

    ! apply the wall boundary condition
    do i=n+1,n+ng
       j = mp_identity(i-n)
       if(fsp(j).eqv..false.)then  !mirror particle only if not mirroring free surface particle
          A(:,j) = A(:,j) + A(:,i)
       end if
       A(:,i) = 0.0
       A(i,:) = 0.0
       B(i) = 0.0
    end do
    
    ! Solve the ppe
    call BiCGStab(A,B,X)  ! no preconditioner

    ! pass X back into P, and put boundary particle values in
    do i=1,n
       if(fsp(i))then
          p(i) = 0.0   ! free surface particle
       else
          p(i) = X(i)  ! regular fluid particle
       end if
    end do
    do i=n+1,n+ng
       p(i) = p(mp_identity(i-n)) ! mirror particle
    end do

    ! calculate dvdt = (-1/ro)*grad(p)
    do i=1,n
       dvdt(i,:) = 0.0
       do k=1,nneigh(i)
          j=neigh_list(i,k)

          ! find some ij-properties
          rij(:) = r(i,:) - r(j,:)
          mrij = sqrt(dot_product(rij,rij))
          dw = dkfunc(mrij,h)
          gradw(:) = dw*rij(:)/max(mrij,1e-6)
          
          dvdt(i,:) = dvdt(i,:) - m(j)*(p(i)+p(j))*gradw(:)/(ro(j)**2)
       end do
    end do          

    deallocate(A)
    deallocate(B)
    deallocate(X)
    
    return
  end subroutine ppe
  
  subroutine BiCGStab(A,b,x)
    ! this subroutine came from the Internet!!
    implicit none
    real,intent(in )                   :: A (:,:)
    real,intent(in )                   :: b ( : )
    real,dimension(1:size(b, dim=1))   :: x
    real,dimension(1:size(b, dim=1))   :: r, rs, v, p, s, t
    real,parameter                     :: e = 1d-5
    real :: rho,rho_prev,alpha,omega,beta,norm_r,norm_b,summesion,temp
    integer :: it=0,err

    if(size(A, dim=1) /= size(A, dim=2)) stop &
         "Error: Improper dimension of matrix A in BiCGStab."

    ! initial guess
    x  = 0.0d0 
    r  = b - matmul(A,x)                           
    rs = r                                          
    rho   = 1.0d0; alpha = 1.0d0; omega = 1.0d0    
    v  = 0.0d0; p  = 0.0d0
    norm_r = sqrt(dot_product(r,r))           
    norm_b = sqrt(dot_product(b,b))                

    do while(norm_r .GT. e*norm_b)                             

       rho_prev = rho                                     
       rho      = dot_product(rs,r)                        
       beta     = (rho/rho_prev) * (alpha/omega)           
       p        = r + beta * (p - omega*v)                 
       v        = matmul(A,p)                              
       alpha    = rho/dot_product(rs,v)                    
       s        = r - alpha*v                              
       t        = matmul(A,s)                              
       omega    = dot_product(t,s)/dot_product(t,t)
       x        = x + alpha*p + omega*s                    
       r        = s - omega*t                              
       norm_r   = norm2(r)!sqrt(dot_product(r,r))
       norm_b   = norm2(b)!sqrt(dot_product(b,b))   
        
       it = it + 1
    end do                                                      
    !print*, "Iteration Required :", it
    !write(6,*) "iterations",it
    
    return
  end subroutine BiCGStab

  subroutine fickian_shift
    real,dimension(n,dims) :: dr
    real,dimension(dims,dims) :: kg_cf
    real,dimension(dims) :: grad_c
    real,dimension(n,dims) :: grad_v1,grad_v2,grad_p
    real :: fij0,fij,rfij,D,dw
    integer :: i,j,k,l
    
    do i=1,n

       ! setup fij
       fij0 = 1.0/kfunc(h/h_ov_spacing,h)

       ! setup kgcf
       call kernel_gradient_correction(i,kg_cf)

       grad_c(:) = 0.0
       grad_v1(:,:) = 0.0
       grad_v2(:,:) = 0.0
       grad_p(:,:) = 0.0

       do k=1,nneigh(i)
          j=neigh_list(i,k)

          ! calculate fij
          fij = kfunc(mrij,h)*fij0
          rfij = 0.2*(abs(min(p(i),0.0))/(ro(i)**2) + abs(min(p(j),0.0))/(ro(j)**2))
          rfij = rfij*fij**4.0

          ! calculate gradw including kgcf
          rij(:) = r(i,:) - r(j,:)
          mrij = sqrt(dot_product(rij,rij))
          dw = dkfunc(mrij,h)
          gradw(:) = dw*rij(:)/max(mrij,1e-6)
          gradw(:) = matmul(kg_cf,gradw)

          ! calculate grad C
          grad_c(:) = grad_c(:) + m(j)*(1.0+rfij)*gradw(:)/ro(j)
          
          ! calculate gradient of velocity.
          grad_v1(i,:) = grad_v1(i,:) - m(j)*(v(i,1)-v(j,1))*gradw(:)/ro(j)
          grad_v2(i,:) = grad_v2(i,:) - m(j)*(v(i,2)-v(j,2))*gradw(:)/ro(j)

          ! calculate grad p
          grad_p(i,:) = grad_p(i,:) - m(j)*(p(i)-p(j))*gradw(:)/ro(j)
          
       end do

       ! calculate D
       D = h*dt*sqrt(dot_product(v(i,:),v(i,:)))
       !D = 0.5*h**2

       ! calculate delta-r
       dr(i,:) = -1.0*D*grad_c(:)
    end do

    ! adjust properties and shift particles...
    do i=1,n
       r(i,:) = r(i,:) + dr(i,:)
       v(i,1) = v(i,1) + dot_product(grad_v1(i,:),dr(i,:))
       v(i,2) = v(i,2) + dot_product(grad_v2(i,:),dr(i,:))
       p(i) = p(i) + dot_product(grad_p(i,:),dr(i,:))       
    end do

    return
  end subroutine fickian_shift

  subroutine calc_ddt_isph
    real :: diss_v(dims)
    integer :: i,j,k,l
    real :: dw

    ! calculate tau if non-newtonian
    if(nnf)then
       call calculate_tau_nnf
    end if
    
    do i=1,n !for all internal particles
       dvdt(i,:) = 0.0 + g(:)
             
       do k=1,nneigh(i) !loop over neighbours
          j=neigh_list(i,k)

          ! calc distance, vdotr and 'ij' properties
          rij(:) = r(i,:) - r(j,:)
          vij(:) = v(i,:) - v(j,:)
          mrij = sqrt(dot_product(rij,rij))
          
          ! find the gradient of kernel function
          dw = dkfunc(mrij,h)
          gradw(:) = dw*rij(:)/max(mrij,1e-6)

          ! calculate the dissipative terms if it's a Non-Newtonian fluid
          if(nnf) then
             diss_v(:) = m(j)*matmul((tau(i,:,:)/ro(i)**2)+(tau(j,:,:)/ro(j)**2),gradw(:))
          else
             diss_v(:) = 2.0*nu*m(i)*dot_product(rij,gradw)*vij(:)/(ro(j)*(mrij**2 + 0.01*h**2))
          end if
          
          !augment forces and
          dvdt(i,:)=dvdt(i,:) + diss_v(:)
       end do
    end do
    
    return
  end subroutine calc_ddt_isph

    subroutine set_dt_isph
    real vsig,minsf,minsc
    integer i,l
    minsc=1e12
    minsf=1e12
    !loop over all internal particles
    do i=1,n
       ! courant constraint
       vsig = dot_product(v(i,:),v(i,:))
       if(h/vsig.le.minsc)then
          minsc = h/vsig
       end if
    end do
    ! viscous constraint
    minsf = h*h/nu

    ! choose the most stringent
    dt=min(0.4*minsc,0.25*minsf)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dt = min(dt,0.01)
    !!!!!!!!!!!!!!!!!!!!!
    return
  end subroutine set_dt_isph

end module evolve_i_mod
