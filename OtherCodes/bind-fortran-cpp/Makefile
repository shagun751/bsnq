
OBJ = testc.o testBind.f90

default:	$(OBJ)
	ifort -o bind $(OBJ) -cpp -lstdc++

testc.o:	testc.cpp
	icc -c testc.cpp