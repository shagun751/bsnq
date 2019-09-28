!!---------------------- Version 1.x.x -----------------------!!
!!  -> Quadratic + Linear
!!  -> Boundary - SemiDirect + Penalty + Gauss Seidel
!!    -> 11 - Inlet - No Absorbing
!!    -> 12 - NoSlip Wall 
!!    -> 13 - Slip Wall - Rectangular wall only
!!    -> 14 - Outlet  - Not coded
!!    -> 15 - Sponge - BOUSS2D approach generalised input
!!  -> Porosity - Energent structure only
!!    -> Generalised input
!!  -> Bottom Shear
!!  -> Turbulence
!!  -> Pressure - v4
!!  -> Solver - Normalised
!!  -> Generalised input - v3
!!  -> GMRES with guess value
!!  -> Drichlet BndCond for eta, p, q
!!  -> Attempting to solve the spurious oscillation issue 
!! 	   in the sponge layers by changing 
!!     the bndType from 12 to 14 (new)
!!  -> Paralution CSR
!!  3.3-> Pressure defined on linear and quadratic nodes
!!	3.5-> XML output
!!	3.6-> Input 3
!!  ! Porosity depth for overtopping (non Sorenson)
!!    -> corrected turb resist coeff to 1.5 from 0.81
!!    -> our model was tuned for 1.5
!!  ! Smoothing Pdt Qdt after some depth to reduce osc
!!
!!	Version 7 ends here
!!
!!  mafi defintions
!!  mafi(1)     Mesh File
!!  mafi(2)     Paraview output
!!  mafi(3)     Volume output
!!  mafi(4)     <Unknown>
!!  mafi(5)     Input file
!!  mafi(6)     Porosity file
!!  mafi(7)     Wave probes files
!!------------------------------------------------------------!!
!!	1.0->
!!---------------------- Version 8.x.x -----------------------!!
!!  Jan 25 2018
!!  1000
!!  continued from bsnq_v7.3.6
!!  1.0   ->  Removed turbulence code
!!        ->  Create code 9 rout file
!!  4.0   -x  Shephard search and cubic spline weight 
!!        ->  Speed calculation using system_clock
!!  8.0   ->  Neumann boundary condition for etadt
!!  21.0  ->  Removing all Kennedy breaking implementation
!!  31.0  ->  Removed smoothing module and some bloatcode

! include 'subroutines/bndIntegral.f90'
! include 'subroutines/boundaryDifferential.f90'
! include 'subroutines/geometry.f90'
! include 'subroutines/gqMatrixSet1.f90'
! include 'subroutines/matrixSet1.f90'
! include 'subroutines/matrixSet2.f90'
! include 'subroutines/mergeSort.f90'
! include 'subroutines/nodeConnAll.f90'
! include 'subroutines/outputN.f90'
! include 'subroutines/porMatrices.f90'
! include 'subroutines/resumeFile.f90'
! include 'subroutines/shipFncs.f90'
! include 'subroutines/waveCalculator_v2.f90'


program boussinesqQuad
use bsnqGlobVars
use bsnqModule
implicit none

interface
  ! Shagun modify 2017_08_14
  subroutine paralution_init(nthreads) BIND(C)
    use, intrinsic :: ISO_C_BINDING, only : C_INT
    integer(kind=C_INT), value, intent(in)  :: nthreads
  end subroutine paralution_init

  subroutine paralution_stop() BIND(C)
  end subroutine paralution_stop    

  subroutine paralution_fortran_solve_coo( n, m, nnz, solver, mformat, preconditioner, pformat,    &
                                          rows, cols, rval, rhs, atol, rtol, div, maxiter, basis, &
                                          p, q, x, iter, resnorm, ierr ) BIND(C)

    use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE, C_CHAR

    integer(kind=C_INT), value, intent(in)  :: n, m, nnz, maxiter, basis, p, q
    real(kind=C_DOUBLE), value, intent(in)  :: atol, rtol, div
    integer(kind=C_INT),        intent(out) :: iter, ierr
    real(kind=C_DOUBLE),        intent(out) :: resnorm
    type(C_PTR),         value, intent(in)  :: rows, cols, rval, rhs
    type(C_PTR),         value              :: x
    character(kind=C_CHAR)                  :: solver, mformat, preconditioner, pformat

  end subroutine paralution_fortran_solve_coo

  subroutine paralution_fortran_solve_csr( n, m, nnz, solver, mformat, preconditioner, pformat, &
                                          ivCSR, jvCSR, rval, rhs, atol, rtol, div, maxiter, basis, &
                                          p, q, x, iter, resnorm, ierr ) BIND(C)

    use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE, C_CHAR

    integer(kind=C_INT), value, intent(in)  :: n, m, nnz, maxiter, basis, p, q
    real(kind=C_DOUBLE), value, intent(in)  :: atol, rtol, div
    integer(kind=C_INT),        intent(out) :: iter, ierr
    real(kind=C_DOUBLE),        intent(out) :: resnorm
    type(C_PTR),         value, intent(in)  :: ivCSR, jvCSR, rval, rhs
    type(C_PTR),         value              :: x
    character(kind=C_CHAR)                  :: solver, mformat, preconditioner, pformat

  end subroutine paralution_fortran_solve_csr
end interface

!!--------------------------Declarations---------------------------!!
  type(bsnqCase)::bq
  type(waveType)::wv
  character(len=C_KSTR)::bqtxt
!!------------------------End Declarations-------------------------!!
  
  call getarg(1,bq%probname)  
  do while(len_trim(bq%probname).lt.1)
    write(*,'(A)',advance='no')"Enter Problem Name:"
    read(*,*)bq%probname
  enddo
  write(*,*)"Problem Name: "//trim(bq%probname)

  bqtxt=trim(bq%probname)//'.rout'
  open(9,file=trim(bqtxt))

  call bq%meshRead
  call bq%femInit
  call bq%setRun  
  call bq%initMat

  call bq%statMatrices  
  call bq%dynaMatrices

  wv=waveType(2d0,10d0,0.5d0,0d0,0d0,90d0)
  write(*,*)wv%T,wv%L,wv%thRad

  write(9,*)"boussinesqQuad End"
end program boussinesqQuad
