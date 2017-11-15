ksph:	./src/particles_mod.f90		\
	./src/output_mod.f90 		\
	./src/neighbours_mod.f90	\
	./src/evolve_i_mod.f90	 	\
	./src/evolve_wc_mod.f90		\
	./src/bound_mod.f90		\
	./src/kernel_mod.f90		\
	./src/viscosity_mod.f90		\
	./src/main.f90
	gfortran -o ksph -g		\
	./src/particles_mod.f90		\
	./src/output_mod.f90 		\
	./src/neighbours_mod.f90	\
	./src/evolve_i_mod.f90	 	\
	./src/evolve_wc_mod.f90		\
	./src/bound_mod.f90		\
	./src/kernel_mod.f90		\
	./src/viscosity_mod.f90		\
	./src/main.f90
clean:
	rm *.mod
clear:
	rm ./output/*.out
