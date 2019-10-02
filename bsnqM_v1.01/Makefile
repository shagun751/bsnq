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

COMPILER=ifort

MAIN=bsnqQuadMain

OBJ=paralution_fortran.o $(MAIN).f90

SUBS=\
subroutines/mods.f90 \
subroutines/bsnqModule.f90 \
subroutines/geometry_v1.0.f90 \
subroutines/mergeSort.f90 \
subroutines/nodeConnAll.f90 \
subroutines/matrixSet1.f90 \
subroutines/matrixSet2.f90 \
subroutines/femAnalyticalTri_v1.0.f90 \
subroutines/bndIntegral.f90 \
subroutines/bndCondition_v1.0.f90 \
subroutines/solver_v1.0.f90
#subroutines/miscOps_v1.0.f90
	
default: $(OBJ)
	mkdir -p Output Resume
	$(COMPILER) -o bsnq $(SUBS) $(OBJ) $(PARA_LIB) $(PARA_FLAGS)
	
# $(MAIN).o: $(MAIN).f90
# 	ifort $(PARA_FLAGS) -c $(MAIN).f90 $(SUBS)

paralution_fortran.o: $(PARA_CONV)
	g++ $(PARA_FLAGS) $(PARA_INC) -c $(PARA_CONV)

clean:
	rm -rf *.dat Output 
	mkdir Output

cleanAll:
	rm -rf *.dat *.o *.mod Output 
	mkdir Output

rmout:
	rm -rf Output
	mkdir Output
	
run: clean default	
	./bsnq rect2D
	cat rect2D.rout


