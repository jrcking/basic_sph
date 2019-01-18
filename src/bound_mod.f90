!! this file is part of Jack King's basic SPH code for simulating a Dam Break
!! 

module bound_mod
  use particles_mod
  implicit none

contains
  subroutine boundaries
    real alpha
    integer nbp
    integer i,l
    integer b_count,b_lh(dims),b_dim(dims)

    ! if a particle has escaped, turn it round and put it back in!
    do i=1,n
       do l=1,dims
          if(r(i,l).lt.rlow(l))then
             r(i,l)= 2.0*rlow(l) - r(i,l)
             v(i,l)=-1.0*v(i,l)
          end if
          if(r(i,l).gt.rhigh(l))then
             r(i,l)= 2.0*rhigh(l) - r(i,l)
             v(i,l)=-1.0*v(i,l)
          end if
       end do
    end do

    ! for every regular particle, find if close to bound
    ! if within 1*h of bound, then will interact with it's ghost
    ! if > 1*h from boundary, won't interact, so no need for ghost.
    alpha = 1.1

    nbp = 0
    do i=1,n
       b_count=0
       do l=1,dims
          ! by a low boundary
          if(r(i,l)-alpha*h.le.rlow(l)) then
             nbp = nbp + 1
             b_count = b_count + 1
             b_lh(b_count) = -1
             b_dim(b_count) = l
             call copy_particle(i,n+nbp)
             call mirror_particle(i,l,n+nbp,-1)
             v(n+nbp,l) = -1.0*v(i,l)
          end if
          ! by a high boundary
          if(r(i,l)+alpha*h.ge.rhigh(l)) then 
             nbp = nbp + 1
             b_count = b_count + 1
             b_lh(b_count) = 1
             b_dim(b_count) = l
             call copy_particle(i,n+nbp)
             call mirror_particle(i,l,n+nbp,1)
             v(n+nbp,l) = -1.0*v(i,l)
          end if
       end do
       ! by a corner - if a particle is "near" 2 or more edges, it
       ! is also near a corner. b_dim and b_lh store info about which
       ! corners it is near.
       if(b_count.eq.2) then
          nbp = nbp + 1
          call copy_particle(i,n+nbp)
          do l=1,b_count
             call mirror_particle(i,b_dim(l),n+nbp,b_lh(l))
             v(n+nbp,b_dim(l)) = -1.0*v(i,b_dim(l))
          end do
       end if
    end do
    ! update ng 
    ng=nbp           

    return
  end subroutine boundaries

  subroutine mirror_particle(i,l,nn,lh)
    integer,intent(in)::i,l,nn,lh

    !adjust position
    if(lh.eq.-1) then
       r(nn,l) = 2.0*rlow(l) - r(i,l)
    else if(lh.eq.1) then
       r(nn,l) = 2.0*rhigh(l) - r(i,l)
    end if    
    return
  end subroutine mirror_particle

  subroutine copy_particle(a,b)
    ! copies properties of particle with index a to
    ! a particle with index b
    integer a,b,l
    v(b,:)=v(a,:)
    r(b,:)=r(a,:)
    p(b)=p(a)
    ro(b)=ro(a)
    c(b)=c(a)
    return
  end subroutine copy_particle
end module bound_mod
