# Develpoment log for bsnqM

## Versions bsnqM_v1.xx.xx

**Details**  
continued from bsnq\_par\_v8.36

- Quad Jacobian = Linear Jacobian for Triangle
- Analytical Integration
- Object Oriented Programming
- Predictor Corrector Time-stepping

-----------------------------------------------

### v1.01

#### Attempting
- Making the code modular so that it can be used by many
- bsnqModule with all the required variables
- It is also required to make it easier to couple with other codes, especially such as MLPG\_R code.
- FEM analytical integrals as general functions so that they can be used later in other code. Check file *femAnalyticalTri_vx.x.f90*

#### List of Work
- [x] Mesh input (Type0)
- [x] FEM Initialisations
- [ ] Wave probes
- [x] Inlet wave characteristics
- [x] Absorbance coeffs
- [ ] Porosity initialisation
- [x] Mass matrices
- [x] Stationary matrix set 1
- [x] Dynamic matrix set 1
- [x] Boundary integrals
- [x] Full momentum mass matrix
- [x] Dirichlet BC
- [x] Conversion to CSR form
- [ ] Neumann BC for eta
- [ ] Time Stepping - Try Adam-bashforth first
- [x] Time Stepping - Predictor Corrector
- [ ] Time Stepping - RK2
- [ ] Bottom shear 
- [ ] Porosity drag terms
- [x] WaveType class with constructor for waveLength
- [x] Paraview XML output

#### Matrices with correct signs
- [x] Mass M1 and M2
- [x] Bs1, Bs2, Bs3, Bs4
- [x] CxF, CyF
- [x] DMat with porosity removed
- [x] Bs5, Bs6
- [x] Advection matrix with porosity removed
- [x] Hydrostatic pressure matrix

#### Modular structure
- bsnqModule
  - type :: bsnqCase
    - procedure ::  initMat
    - procedure ::  meshRead
    - procedure ::  femInit
    - procedure ::  setRun
    - procedure ::  statMatrices
    - procedure ::  dynaMatrices
    - procedure ::  destructR1
  - type :: waveType
    - constructor :: waveLenCalc
- bsnqGlobVars  
  - Datatypes and constants only

-----------------------------------------------

### References
1. Sørensen, O. R., Schäffer, H. A., & Sørensen, L. S. (2004). Boussinesq-type modelling using an unstructured finite element technique. Coastal Engineering, 50(4), 181–198. [DOI](https://doi.org/10.1016/j.coastaleng.2003.10.005)

1. Agarwal, S., Sriram, V., & Murali, K. (2019). Modelling Wave Interaction with Porous Structures Using Boussinesq Equations. In Proceedings of the Fourth International Conference in Ocean Engineering (ICOE2018) (pp. 573–583). <https://doi.org/10.1007/978-981-13-3119-0_35>

