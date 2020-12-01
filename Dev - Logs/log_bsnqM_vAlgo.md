## Common algorithm developments in the code

1. [Paralution solver FORTRAN plugin [2020-10-15]](#log_bsnqM_vAlgo_1)


### Attempting
- From 2020-Dec-01 onwards, the common algorithm developments will be noted in this document.
- This log file should only be updated for a branch if the bug-fix has been done in that branch.


### List of Work


-----------------------------------------------


<a name = 'log_bsnqM_vAlgo_1' ></a>

### Paralution solver FORTRAN plugin [2020-10-15]
- The existing system of using the paralution solver in FORTRAN was provided by them itself and was robust and simple.
- Paralution was written in C++
- Our code is in FROTRAN
- The coupling between FORTRAN and C++ was limited to using common datatypes.
- The way the existing code was written, every-time we called _paralution_fortran_solve_csr()_ subroutine in FORTRAN, it would do the following
	1. Create a copy of all the A, X, B in C++
	2. Initialise a C++ solver object 
		1. Solver type
		2. Solver tolerance limits and other setting
		3. Save the large A matrix into the solver object (Very slow process)
		4. Save the B and X0 in the solver object
	3. Solve using the settings provided by the user
	4. Copy the solved X from C++ object to the FORTRAN memory
	5. Delete the C++ object.
- The above process happens at each call, which is 12 times per time-step
- In Bsnq the A matrix never changes. Also the solver setting never change. 
- Only B and X0 change
- So the time taken in Steps 1-2 is useless.
- Ideally we would like to create the solver object one time at the beginning of the simulation and save the A matrix into the solver. After that only update B and run the solver to get the solution.
- However the big problem is that the solver is a C++ object and hence cannot be "saved" as a FORTRAN object.
- In the existing code the scope of the solver object is limited to _paralution_fortran_solve_csr()_ function, and as soon as it finishes execution the object space will be deleted and will have to be re-initialised.
- **So the challenge was to figure out how to create and save a C++ object from the FORTRAN code**
- This is where the magic of **void** pointer in C++ comes into picture. **This is the greatest feature of a programming language I have ever seen**. Refer to this [link](https://community.intel.com/t5/Intel-Fortran-Compiler/Calling-C-cpp-objects-from-a-Fortran-subroutine/td-p/1110556)

#### Modifications to the solver code
- I have moved the plugin code from the general location (where paralution is installed) to the subroutines folder instead.
- The file name is still the same _paralution_fortran.cpp_. This way one can switch to the old generic plugin just by changing the makefile
- The modification are quite intensive and difficult to follow.
- To understand the implementation, first refer to the code in 'OtherCode/bind-fortran-cpp'. Looking directly at the modifications in the Bsnq code may be very confusing
- This is how the process works for each.
	1. Create the desired class in C++.
	2. Create a _void_ pointer in C++.
	3. Create a function in C++ which will be called from FORTRAN. 
		- This function will create the solver object in C++ and return a void pointer pointing to this object. 
		- As pointer type is common between C++ and FORTRAN, as long as this pointer is alive the object that it is pointing to will be alive.
		- This way we can create a C++ object and point to it using FORTRAN pointer
		- The interesting thing to note is that pointers are generally of the datatype at which they are pointing
		- However void pointer is special as it can point to any datatype and hence it can point to object of any class. (beautiful feature)
	4. Now create function in C++ which can be called from FORTRAN. 
		- Pass the void pointer from FORTRAN in these functions to refer to a specific object along with the other variables. 
		- And then use the normal C++ code to call the function associate with that object with the respective variable in C++ itself.
- Its a very complex process but it works well

The new solution process in bsnq code is as below
	- The process is now modified where a solver object is create for W, Eta and PQ equations each.
	- During the creation the solver setting are fixed to BiCGStab
	- Error norm is fixed to L2 norm
	- All tolerences are set once in the beginning
	- The large A matrix is transferred to the solver only once in the beginning and after that it is never transferred again. This is where I think we gain the majority of the time in transferring A and allocating space for it. In old method this method happened at each call.
	- Space is create for X and B
	- Only X0 and B are updated every time-step and solver is just run (**not re-created every time**)


In the Boussinesq code modification are made to

1. subroutines/bsnqModule.f90
	- Creating the solver object for W, Eta and PQ equations separately
2. subroutines/bsnqModuleFncs.f90
	- Changed the solver calling lines. Old lines are commented to easily switch back to old reliable robust version if required.
3. subroutines/solver_v1.0.f90
	- New solver subroutine. _solveSys2()_
	- Old one is also there but is un-used _solveSys()_
4. subroutines/paralution_fortran.cpp
	- The new file with my custom plugin code.

Although I think this is the best piece of code I have ever written, it gave very limited gains in speed.

#### Speed gains
- Overall on normal i7-9th Gen system I got upto 1.08x
- However on Aqua server it hardly made any significant difference.

-----------------------------------------------