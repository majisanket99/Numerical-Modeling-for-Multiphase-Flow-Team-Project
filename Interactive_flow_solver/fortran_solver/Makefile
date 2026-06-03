# Fortran compiler
FC = gfortran

# Compiler flags
FFLAGS = -O3 -w -fallow-argument-mismatch

# Source files
FILES = module_globalvariables.F90      \
        module_fluid.F90                \
        module_immersedboundary.F90     \
        module_floodfill.F90            \
        initialize.F90 output.F90       \
        solve.F90                       \
        deallocation.F90                \
        IBM2D.F90

# Build FFT library
fft:
	$(FC) -c fftpack5.1.f90 $(FFLAGS)

# Build IBM solver
ibm2d:
	$(FC) fftpack5.1.o $(FILES) $(FFLAGS) -o ibm2d.out

all: fft ibm2d

# Clean generated files
clean:
	rm -f *.o *.mod *.MOD ibm2d.out
	rm -f ./results/*
	rm -f xstop
