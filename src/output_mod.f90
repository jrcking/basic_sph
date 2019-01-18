!! this file is part of Jack King's basic SPH code for simulating a Dam Break
!! 

module output_mod
  use particles_mod
  implicit none
contains

  subroutine create_files
    open(unit=21,file="./output/dt.out")
    open(unit=22,file="./output/nbp.out")
    return
  end subroutine create_files

  subroutine output_simulation_data(n)
    integer,intent(in) :: n
    write(21,*) n,dt
    write(22,*) n,ng
    return
  end subroutine output_simulation_data

  subroutine close_files
    close(21)
    close(22)
    return
  end subroutine close_files
  
  subroutine output_particle_data(n_out)
    integer :: i,n_out
    character :: tmp*(5)

    ! put n_out in a character variable
    write(tmp,'(i5.5)') n_out
    tmp=adjustl(tmp)

    ! open the files
    open(unit=9,file="./output/r"//tmp//".out")
    open(unit=10,file="./output/v"//tmp//".out")
    open(unit=11,file="./output/p"//tmp//".out")
    open(unit=12,file="./output/ro"//tmp//".out")
    open(unit=13,file="./output/neighbours"//tmp//".out")

    ! write position and velocity files
    ! assume dims=2 (because I have no intention of going 3d)
    do i=1,n!+ng
       write(9,*) r(i,1),r(i,2)
       write(10,*) v(i,1),v(i,2)
    end do
 
    ! write  density, pressure and # neighbours to files
    do i=1,n!+ng
       write(11,*) p(i)
       write(12,*) ro(i)
       write(13,*) nneigh(i)
    end do
    
    !close the files
    do i=9,13
       close(i)
    end do
    
    return
  end subroutine output_particle_data

end module output_mod
