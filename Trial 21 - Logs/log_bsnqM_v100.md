# Develpoment log for bsnqM

## Versions bsnqM_v1.xx.xx

**Details** \
continued from bsnq_par_v8.36
- Quad Jacobian = Linear Jacobian for Triangle
- Analytical Integration
- Object Oriented Programming
- Predictor Corrector Time-stepping
-----------------------------------------------

### v1.01

##### Attempting
- bsnqModule with all the required variables

##### List of Work
- [ ] Wave probes
- [ ] Inlet wave characteristics
- [ ] Absorbance coeffs
- [ ] Porosity initialisation
- [x] Stationary matrix set 1
- [ ] Stationary matrix set 2

##### Matrices with correct signs
- [x] Bs1, Bs2, Bs3, Bs4
- [x] CxF, CyF
- [x] DMat with porosity removed
- [ ] Bs5, Bs6

-----------------------------------------------

### References
[1] Sørensen, O. R., Schäffer, H. A., & Sørensen, L. S. (2004). Boussinesq-type modelling using an unstructured finite element technique. Coastal Engineering, 50(4), 181–198. https://doi.org/10.1016/j.coastaleng.2003.10.005

[2] Agarwal, S., Sriram, V., & Murali, K. (2019). Modelling Wave Interaction with Porous Structures Using Boussinesq Equations. In Proceedings of the Fourth International Conference in Ocean Engineering (ICOE2018) (pp. 573–583). https://doi.org/10.1007/978-981-13-3119-0_35

