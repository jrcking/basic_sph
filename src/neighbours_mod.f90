!! this file is part of Jack King's basic SPH code for simulating a Dam Break
!! 

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

end module neighbours_mod
