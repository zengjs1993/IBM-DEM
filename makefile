#f90 =/usr/local/mpi/bin/mpif90
#f90=ifort
f90=mpif90


#mpif90 -c -fPIC -Wall -ffree-line-length-0 -Wno-unused-dummy-argument -g    -I/gpfs/share/home/1301111457/petsc-3.10.2/include -I/gpfs/share/home/1301111457/petsc-3.10.2/arch-linux2-c-debug/include    -o ex45f3d.o ex45f3d.F90
#mpif90 -fPIC -Wall -ffree-line-length-0 -Wno-unused-dummy-argument -g   -o ex45f3d ex45f3d.o  -Wl,-rpath,/gpfs/share/home/1301111457/petsc-3.10.2/arch-linux2-c-debug/lib -L/gpfs/share/home/1301111457/petsc-3.10.2/arch-linux2-c-debug/lib -Wl,-rpath,/gpfs/share/home/1301111457/petsc-3.10.2/arch-linux2-c-debug/lib -L/gpfs/share/home/1301111457/petsc-3.10.2/arch-linux2-c-debug/lib -Wl,-rpath,/gpfs/share/software/compilers/openmpi/3.0.0/gcc/4.8.5/lib -L/gpfs/share/software/compilers/openmpi/3.0.0/gcc/4.8.5/lib -Wl,-rpath,/usr/lib/gcc/x86_64-redhat-linux/4.8.5 -L/usr/lib/gcc/x86_64-redhat-linux/4.8.5 -lpetsc -lflapack -lfblas -lm -lX11 -lstdc++ -ldl -lmpi_usempi -lmpi_mpifh -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lpthread -lstdc++ -ldl
#PETSC_DIR=/home/junsheng/petsc-3.10.2
#include ${PETSC_DIR}/lib/petsc/conf/variables
#include ${PETSC_DIR}/lib/petsc/conf/petscvariables
#include ${PETSC_DIR}/lib/petsc/conf/rules
#opt = -fast
opt=-O3

srs=vars.F90 ZJS_VARS.F90 main.F90 boundary.F90 io.F90 initialization.F90 allocate.F90 poisson.F90 fluid.F90 interpolation.F90 ZJS_AUXILIARY_DEM.F90 ZJS_COLLISION_DEM.F90 ZJS_DRAGFORCE_DEM.F90 ZJS_EVOLUTION_DEM.F90 ZJS_INIT_DEM.F90 ZJS_INTERPOLATE_DEM.F90 ZJS_IO_DEM.F90 ZJS_MPI_DEM.F90 ZJS_UPDATE_LIST_DEM.F90 ZJS_VOLUMEFRACTION_DEM.F90 ZJS_MS_DEM.F90 couplesolverCylinder.F90

OBJS=$(srs:.F90=.o)

%.o:%.F90
	$(f90) $(opt) -c -fPIC -Wall -ffree-line-length-0 -Wno-unused-dummy-argument $<   -I/gpfs/share/home/1301111457/petsc-3.10.2/include -I/gpfs/share/home/1301111457/petsc-3.10.2/arch-linux2-c-debug/include 
default: $(OBJS)
	$(f90) -O3 -o ibm-dem $(OBJS) -g -Wl,-rpath,/gpfs/share/home/1301111457/petsc-3.10.2/arch-linux2-c-debug/lib -L/gpfs/share/home/1301111457/petsc-3.10.2/arch-linux2-c-debug/lib -Wl,-rpath,/gpfs/share/home/1301111457/petsc-3.10.2/arch-linux2-c-debug/lib -L/gpfs/share/home/1301111457/petsc-3.10.2/arch-linux2-c-debug/lib -Wl,-rpath,/gpfs/share/software/compilers/openmpi/3.0.0/gcc/4.8.5/lib -L/gpfs/share/software/compilers/openmpi/3.0.0/gcc/4.8.5/lib -Wl,-rpath,/usr/lib/gcc/x86_64-redhat-linux/4.8.5 -L/usr/lib/gcc/x86_64-redhat-linux/4.8.5 -lpetsc -lflapack -lfblas -lm -lX11 -lstdc++ -ldl -lmpi_usempi -lmpi_mpifh -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lpthread -lstdc++ -ldl
clean:
	rm -f ibm-dem *.mod *.o
