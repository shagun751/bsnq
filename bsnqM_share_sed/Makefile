IFORTCHK=-check all -fpe0 -warn -traceback

COMPILER=ifort
#COMPILER=gfortran

TIMESTEPPING = -DTRK4      #RK4 TimeStepping
#TIMESTEPPING = -DTADBAE3  #Adam-Bashforth Explicit 3 point

ifeq ($(COMPILER),ifort)	
	C_COMP=icc
	PARA_FLAGS=-O3 -lstdc++ -fopenmp -L/usr/local/cuda/lib64	
else ifeq ($(COMPILER),gfortran)	
	C_COMP=g++
	PARA_FLAGS=-O3 -lstdc++ -fopenmp -L/usr/local/cuda/lib64	
endif
OTH_FLAGS = -cpp $(TIMESTEPPING) -I./lib/

PARA_LIB=./lib/libparalution.so 
PARA_INC=-I$(PARA_PATH)/build/inc

MAIN=bsnqQuadMain
MODULE=bsnqModule

OBJ= lib/paralution_fortran.o lib/$(MODULE).o $(MAIN).f90
	
default: $(OBJ)
	mkdir -p Output
	$(COMPILER) -o bsnq $(OBJ) $(PARA_LIB) $(PARA_FLAGS) $(OTH_FLAGS)
	
clean:
	rm -rf Output 
	mkdir Output

cleanAll:
	rm -rf *.o Output bsnq 
	mkdir Output

rmout:
	rm -rf Output
	mkdir Output
	
run: clean default	
	./bsnq rect2D
	cat rect2D.rout


