# Develpoment log for bsnqM

## Versions bsnqM v1.01

**Details**  
continued from bsnq\_par\_v8.36

- Quad Jacobian = Linear Jacobian for Triangle
- Analytical Integration
- Object Oriented Programming
- Modular Time-stepping formulation
  
1. [Initial Development log](./log_bsnqM_v0001.md)
1. [Moving pressure field and Gradient MLS development](./log_bsnqM_v0002.md)
1. [Vertical velocity calculation](#log_bsnqM_v0003)

-----------------------------------------------

<a name = 'log_bsnqM_v0003' />

## Vertical velocity calculation

### Attempting
- Calculate velocities along the depth 


### List of Work
- [x] Consecutive derivative based - uDx uDxx uDxxx 
- [x] Consecutive derivative based - pDx pDxx pDxxx 


### Observations : 3rd Derv consecutive derivative based [2020-03-06]
- Analytical function sin(x) was used to check till third derivative calculated using the MLS code.
- The consecutive derivative based approach is:
	- f'(x) = d ( f(x) )/ dx
	- f''(x) = d ( f'(x) )/ dx
	- f'''(x) = d ( f''(x) )/ dx  
- The subroutine is _calcDerv_, part of bsnqModule. It will only be called if the pObf is allocated by the subroutine _setMFree_, also a part of bsnqModule.
- From the results below it can be seen that near the boundaries the 2nd and 3rd derivatives are inaccurate. **This can be a issue in calculating depth resolved velociies for coupling near the boundaries.**

<p align="centre">  <img width="90%" src="./mlsVsAnalitcalDervSinx/Derv_1.jpg">

<img width="90%" src="./mlsVsAnalitcalDervSinx/Derv_2.jpg">

<img width="90%" src="./mlsVsAnalitcalDervSinx/Derv_3.jpg">  

**Fig :** Results of 1st, 2nd and 3rd derivatives of sin(x) compared for MLS code against.
</p>

-----------------------------------------------

## References
1. Bosboom, Judith. 1995. “Boussinesq Modelling of Wave-Induced Particle Velocities.” TUDelft.

1. Dingemans, M. W. 1994. “Water Wave Propagation over Uneven Bottoms.” TUDelft.
