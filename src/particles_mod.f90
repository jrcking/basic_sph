!! this file is part of Jack King's basic SPH code for simulating a Dam Break
!! 

module particles_mod
  use kernel_mod
  implicit none
  integer,parameter :: dims = 2
  real, allocatable,dimension(:) :: ro,p,c,drodt
  real, allocatable,dimension(:,:) :: v,r,dvdt
  integer,allocatable,dimension(:) :: nneigh
  integer,allocatable,dimension(:,:) :: neigh_list
  real,dimension(dims) :: rlow,rhigh,g,pb_low,pb_high
  integer,dimension(dims) :: nx
  real,dimension(dims) :: rij,vij,gradw
  real :: mrij,p_space
  real :: dt,tmax,dmp_period
  real :: m,h,ro_init,av_value,h_ov_spacing
  integer :: n,ng,nb
  logical :: slip

contains
  subroutine create_domain
    integer :: l

    ! read in run parameters
    open(unit=1,file="simulation_parameters",status='old')
    read(1,*) nx
    n = nx(1)*nx(2)
    ng = n
    read(1,*) tmax
    read(1,*) dmp_period
    read(1,*) av_value
    read(1,*) ro_init
    read(1,*) rlow
    read(1,*) rhigh
    read(1,*) g
    read(1,*) h_ov_spacing
    read(1,*) pb_low
    read(1,*) pb_high
    close(1)
 
    return
  end subroutine create_domain
  
  subroutine create_particles
    real,dimension(dims)::dx,rij,tmp_x
    real :: mrij,p_space,tmp_circle
    integer,dimension(dims) :: tmp2(dims)
    integer :: i,j,l,nc
    
    ! create arrays for particle properties
    allocate(ro(n+ng))
    allocate(v(n+ng,dims))
    allocate(r(n+ng,dims))
    allocate(drodt(n))
    allocate(dvdt(n,dims))
    allocate(p(n+ng))
    allocate(c(n+ng))
    allocate(nneigh(n+ng))
    allocate(neigh_list(n+ng,n+ng))

    ! # particles in each direction
    dx(:)=(pb_high(:)-pb_low(:))/real(nx(:))
    p_space = dx(1)

    ! place the particles
    do i=1,n
       j=floor(real(i-1)/real(nx(1)))
       tmp2(1)=mod(i-1,nx(1))
       tmp2(2)=mod(j,nx(2))
       r(i,:) = pb_low(:) + 0.5*dx(:) + dx(:)*tmp2(:)
    end do

    ! set other particle properties
    ro(:) = ro_init
    v(:,:) = 0.0
    ! v(:,1) = r(:,2)*(2.0*pb_high(2)-r(:,2))*64.0
    h = h_ov_spacing*p_space
    m = ro_init*(h/h_ov_spacing)**dims
    
    ! write numbers of particles
    write(6,*) n,ng
    
    return
  end subroutine create_particles
  
  subroutine destroy_particles
    deallocate(ro)
    deallocate(v)
    deallocate(r)
    deallocate(drodt)
    deallocate(dvdt)
    deallocate(p)
    deallocate(c)
    deallocate(nneigh)
    deallocate(neigh_list)
    return
  end subroutine destroy_particles

end module particles_mod
