!!---------------------- Version 1.x.x -----------------------!!
!!  -> Quadratic + Linear
!!  -> Boundary - SemiDirect + Penalty + Gauss Seidel
!!    -> 11 - Inlet - No Absorbing
!!    -> 12 - NoSlip Wall 
!!    -> 13 - Slip Wall - Rectangular wall only
!!    -> 14 - Sponge - BOUSS2D approach generalised input
!!    -> 15 - Outlet  - Not coded
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


program boussinesqQuad
use bsnqGlobVars
use bsnqModule
implicit none

!!--------------------------Declarations---------------------------!!
  real(kind=C_K2)::rTime

  type(bsnqCase)::bq  
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

  call system_clock(bq%sysC(1))
  call bq%meshRead
  call bq%femInit
  call bq%setRun  
  call bq%initMat
  call bq%statMatrices  
  
  do while(bq%tOb(0)%rtm .lt. bq%endTime)  

    call bq%preInstructs

    !!-------------RK2 S1--------------!!
    bq%ur = bq%tOb(1)%p / bq%tOb(1)%tD
    bq%vr = bq%tOb(1)%q / bq%tOb(1)%tD
    bq%pbpr = bq%tOb(1)%p / bq%por
    bq%qbpr = bq%tOb(1)%q / bq%por   
    rTime=bq%tOb(1)%rtm
    call bq%dynaMatrices(bq%tOb(1)%tD, bq%ur, bq%vr)
    call bq%solveAll(rTime, bq%tOb(1)%p, bq%tOb(1)%q, &
      bq%pbpr, bq%qbpr, bq%tOb(1)%e, &
      bq%gXW, bq%gXE, bq%gXPQ, bq%gRE, bq%gRPQ)
    call bq%updateSoln(1)
    !!-----------End RK2 S1------------!!

    !!-------------RK2 S2--------------!!
    bq%ur = bq%tOb(0)%p / bq%tOb(0)%tD
    bq%vr = bq%tOb(0)%q / bq%tOb(0)%tD
    bq%pbpr = bq%tOb(0)%p / bq%por
    bq%qbpr = bq%tOb(0)%q / bq%por 
    rTime=(bq%tOb(0)%rtm + bq%tOb(1)%rtm)/2d0
    call bq%dynaMatrices(bq%tOb(0)%tD, bq%ur, bq%vr)
    call bq%solveAll(rTime, bq%tOb(0)%p, bq%tOb(0)%q, &
      bq%pbpr, bq%qbpr, bq%tOb(0)%e, &
      bq%gXW, bq%gXE, bq%gXPQ, bq%gRE, bq%gRPQ)
    call bq%updateSoln(2)
    !!-----------End RK2 S2------------!!

    call bq%postInstructs

  enddo

  call system_clock(bq%sysC(2))
  write(9,*)"[MSG] boussinesqQuad End"
  write(9,'(" [TIM] ",F15.4)')1d0*(bq%sysC(2)-bq%sysC(1))/bq%sysRate
  close(9)
end program boussinesqQuad