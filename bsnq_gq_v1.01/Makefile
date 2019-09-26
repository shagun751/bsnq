SWITCH_ICC_RELEASE=-ipo -O3 -vec_report1
SWITCH_ICC_CUSTOM=-O3
SWITCH_ICC_DEBUG= -g -debug all -check all -implicitnone -fp-stack-check -heap-arrays -ftrapuv -check pointers -check bounds -traceback
SWITCH_GFORTRAN_RELEASE=-Ofast
SWITCH_GFORTRAN_O2=-O2
SWITCH_GFORTRAN_DEBUG=-O0 -g -fbounds-check -fbacktrace

PARA_PATH=$(HOME)/Programs/paralution_gcc_min
PARA_FLAGS=-O3 -lstdc++ -fopenmp -L/usr/local/cuda/lib64
PARA_CONV=$(PARA_PATH)/src/plug-ins/paralution_fortran.cpp
PARA_LIB=$(PARA_PATH)/build/lib/libparalution.so
PARA_INC=-I$(PARA_PATH)/build/inc

MAIN=bsnqQuadMain

OBJ=paralution_fortran.o $(MAIN).o
	
default: $(OBJ)
	mkdir -p Output Resume
	gfortran -o bsnq $(OBJ) $(PARA_LIB) $(PARA_FLAGS) -Wall

$(MAIN).o: $(MAIN).f90
	gfortran $(PARA_FLAGS) -c $(MAIN).f90

mmio.o: mmio.f
	gfortran $(PARA_FLAGS) -c mmio.f

paralution_fortran.o: $(PARA_CONV)
	g++ $(PARA_FLAGS) $(PARA_INC) -c $(PARA_CONV)

clean:
	rm -rf *.log *.var *.vtk *.dat *.o Output vol.png fort.20
	mkdir Output

rmout:
	rm -rf Output
	mkdir Output
	



