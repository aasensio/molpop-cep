rem Cleanup:
del *.o
del *.mod

@echo Done with cleanup; now to compilation
@pause

gfortran -c global_molpop.f90
gfortran -c maths_molpop.f90 
gfortran -c coll_molpop.f90
gfortran -c cep_molpop_interface.f90
gfortran -c io_molpop.f90
gfortran -c maths_cep.f90	
gfortran -c constants_cep.f90
gfortran -c global_cep.f90	
gfortran -c io_cep.f90
gfortran -c functions_cep.f90
gfortran -c -DNO -x f95-cpp-input escape_cep.f90
gfortran -c sol_cep.f90
gfortran -c sol_molpop.f90
gfortran -c molpop.f90 

gfortran coll_molpop.o global_molpop.o io_molpop.o maths_molpop.o molpop.o sol_molpop.o cep_molpop_interface.o sol_cep.o constants_cep.o	global_cep.o maths_cep.o escape_cep.o functions_cep.o io_cep.o -o molpop.exe

@echo Done!
@pause
