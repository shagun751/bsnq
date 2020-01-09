COMPILER=gfortran

MAIN=fnpt2Bsnq
MODULE=fnptModule

OBJ=$(MODULE).o $(MAIN).f90

SUBS=\
$(MODULE).for \

default: $(OBJ)
	$(COMPILER) -o fnpt2Bsnq $(OBJ) 
	
$(MODULE).o: $(SUBS)
	$(COMPILER) -c $<

clean:
	rm -rf *.dat Output 

cleanAll:
	rm -rf *.dat *.o subroutines/*.mod *.mod Output bsnq 



