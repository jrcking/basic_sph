module viscosity_mod
  use particles_mod
  use neighbours_mod
  use kernel_mod
  implicit none

contains

    subroutine calculate_tau_nnf
    ! As in Xenakis et al 2015 and Xenakis et al 2015
    ! but I use the kernel to calculate the velocity gradients, rather than FD
    real, dimension(dims) :: gradu,gradv
    real, dimension(dims,dims) :: D,kg_cf
    real :: mu_eff,D_spi
    integer :: i,j,k,l

    ! calculate tau
    do i=1,n+ng

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
       
       gradu(:) = 0.0
       gradv(:) = 0.0
       
       do k=1,nneigh(i)
          j=neigh_list(i,k)
          rij(:) = r(i,:) - r(j,:)
          vij(:) = v(i,:) - v(j,:)
          mrij = sqrt(dot_product(rij,rij))
          gradw(:) = dkfunc(mrij,h)*rij(:)/max(mrij,1e-6)
          gradw(:) = matmul(kg_cf,gradw)

          gradu(:) = gradu(:) - m(j)*vij(1)*gradw(:)/ro(j)
          gradv(:) = gradv(:) - m(j)*vij(2)*gradw(:)/ro(j)          
       end do
       D(1,1) = 2.0*gradu(1)
       D(1,2) = gradu(2) + gradv(1)
       D(2,1) = D(1,2)
       D(2,2) = 2.0*gradv(2)
       D_spi = sqrt(0.5*(D(1,1)**2 + 2*D(1,2)**2 + D(2,2)**2))
       
       ! Rheological model:
       !call rheology_newtonian(D_spi,mu,mu_eff)
       call rheology_power_law(D_spi,mu,mu_eff)
       !call rheology_bingham(D_spi,mu,mu_eff)
       
       tau(i,:,:) = mu_eff*D(:,:)
    end do
  
    return
  end subroutine calculate_tau_nnf

  subroutine rheology_newtonian(SR,mu,mu_eff)
    real,intent(in)::SR,mu
    real,intent(out)::mu_eff

    mu_eff = mu
  end subroutine rheology_newtonian

  subroutine rheology_power_law(SR,mu,mu_eff)
    real,intent(in)::SR,mu
    real,intent(out)::mu_eff
    integer :: n

    ! power law exponent: <1=shear thinning; >1=shear thickening
    n = 0.15 !0.15
    
    mu_eff = mu*max(SR,1e-6)**(n - 1.0)  ! power law
  end subroutine rheology_power_law
  
  subroutine rheology_bingham(SR,mu,mu_eff)
    real,intent(in)::SR,mu
    real,intent(out)::mu_eff
    real :: tau_yield,alpha
    
    ! parameters of bingham fluid
    tau_yield = 25.0
    alpha = 1e2

    if (SR.le.tau_yield/(mu*alpha)) then
       mu_eff = mu*alpha
    else
       mu_eff = mu + tau_yield/SR
    end if
  end subroutine rheology_bingham

end module viscosity_mod
