module neighbours_mod
  use particles_mod
  use kernel_mod
  implicit none

contains  
  subroutine find_neighbours
    real :: n_rad
    integer :: i,j

    ! this is not the most efficient way of doing this...
    n_rad=2.0
    do i=1,n+ng
       nneigh(i)=0
       do j=1,n+ng ! look at every other particle
          rij(:) = r(i,:) - r(j,:)
          mrij=sqrt(dot_product(rij,rij))
          !test for proximity, and add to neighbour list
          if(mrij.le.n_rad*h)then
             nneigh(i)=nneigh(i)+1
             neigh_list(i,nneigh(i))=j
          end if
       end do
    end do
    return
  end subroutine find_neighbours

  subroutine kernel_gradient_correction(i,kg_cf)
    ! As in Gomez-Gesteira et al, 2012, section 2.8.2
    integer,intent(in) :: i
    real,intent(out),dimension(dims,dims) :: kg_cf
    real :: M11,M12,M22,detM,dw
    integer :: j,k

    M11 = 0.0
    M12 = 0.0
    M22 = 0.0
    do k=1,nneigh(i)
       j=neigh_list(i,k)
       rij(:) = r(i,:) - r(j,:)
       mrij = max(sqrt(dot_product(rij,rij)),1e-5) !avoid singularity when i=j
       dw = dkfunc(mrij,h)
       ! elements of M (only 3 as it is symmetric
       M11 = M11 - m(j)*dw*rij(1)*rij(1)/(ro(j)*mrij)
       M12 = M12 - m(j)*dw*rij(1)*rij(2)/(ro(j)*mrij)
       M22 = M22 - m(j)*dw*rij(2)*rij(2)/(ro(j)*mrij)       
    end do
    ! find the inverse of M and allocate into kg_cf
    detM = 1/(M11*M22 - M12**2.0)
    kg_cf(1,1) = M22/detM
    kg_cf(1,2) = -1.0*M12/detM
    kg_cf(2,1) = kg_cf(1,2)-1.0*M12/detM
    kg_cf(2,2) = M11/detM
    return
  end subroutine kernel_gradient_correction

end module neighbours_mod
