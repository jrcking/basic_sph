program main
  use particles_mod
  use output_mod
  use neighbours_mod
  use evolve_wc_mod
  use evolve_i_mod
  use bound_mod
  use kernel_mod
  implicit none
  real :: t
  integer :: i,nt,nt_out

  ! create particles
  call create_domain
  call create_particles

  ! initial setup of boundaries, neighbours and pressure, d/dt
  call calc_p_from_ro
  call boundaries
  call find_neighbours
  call calc_ddt_wcsph

  !some final initialisation
  call create_files
  t=0.0
  nt = 0
  nt_out = 0

  ! the time loop
  do while(t.le.tmax)
     nt = nt + 1
!     do nt=1,200

     ! write the output every dmp_period
     if(floor(t/dmp_period).ge.nt_out) then
        nt_out = nt_out+1
        call output_particle_data(nt_out)
     end if
     call output_simulation_data(nt)
     
     ! evolve the particles
     if(compressible)then
        call evolve_wcsph(nt)
     else
        call evolve_isph(nt)
     end if
     
     ! write some things to screen
     t=t+dt
     write(6,*) nt,t,dt
  end do

  ! end of run housekeeping
  call close_files
  call destroy_particles

  return
end program main
