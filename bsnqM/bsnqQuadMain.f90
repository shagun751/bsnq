!!------------------------ Version M 1.x.x ------------------------!!
!!  -> Quadratic + Linear
!!  -> Boundary - SemiDirect + Penalty + Gauss Seidel
!!    -> 11 - Inlet - No Absorbing
!!    -> 12 - NoSlip Wall 
!!    -> 13 - Slip Wall - Rectangular wall only
!!    -> 14 - Sponge - BOUSS2D approach generalised input
!!    -> 15 - Outlet  - Not coded
!!  -> Porosity - Emergent structure only (not done yet)
!!    -> Generalised input
!!  -> Solver - Normalised
!!  -> Generalised input - v3
!!  -> GMRES
!!  -> Dirichlet BndCond for eta, p, q
!!  -> Paralution CSR
!!	-> XML output
!!-----------------------------------------------------------------!!
!! Time-Stepping : RK4 

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
  call bq%setMFree 
  call bq%initMat  
  call bq%statMatrices    


  do while(abs(bq%tOb(0)%rtm-bq%endTime).gt.bq%dt/2d0)
  
    call bq%preInstructs

    !!-------------RK4 S1--------------!!
    bq%ur = bq%tOb(1)%p / bq%tOb(1)%tD
    bq%vr = bq%tOb(1)%q / bq%tOb(1)%tD
    bq%pbpr = bq%tOb(1)%p / bq%por
    bq%qbpr = bq%tOb(1)%q / bq%por   
    rTime=bq%tOb(1)%rtm
    call bq%dynaMatrices(rTime,bq%tOb(1)%tD, bq%ur, bq%vr)
    call bq%solveAll(rTime, bq%tOb(1)%p, bq%tOb(1)%q, &
      bq%pbpr, bq%qbpr, bq%presr, bq%tOb(1)%e, &
      bq%gXW, bq%gXE, bq%gXPQ, bq%gRE, bq%gRPQ, bq%sysC)
    call bq%updateSoln(1)
    !!-----------End RK4 S1------------!!

    !!-------------RK4 S2--------------!!
    bq%ur = bq%tOb(0)%p / bq%tOb(0)%tD
    bq%vr = bq%tOb(0)%q / bq%tOb(0)%tD
    bq%pbpr = bq%tOb(0)%p / bq%por
    bq%qbpr = bq%tOb(0)%q / bq%por   
    rTime=(bq%tOb(0)%rtm + bq%tOb(1)%rtm)/2d0
    call bq%dynaMatrices(rTime,bq%tOb(0)%tD, bq%ur, bq%vr)
    call bq%solveAll(rTime, bq%tOb(0)%p, bq%tOb(0)%q, &
      bq%pbpr, bq%qbpr, bq%presr, bq%tOb(0)%e, &
      bq%gXW, bq%gXE, bq%gXPQ, bq%gRE, bq%gRPQ, bq%sysC)
    call bq%updateSoln(2)
    !!-----------End RK4 S2------------!!

    !!-------------RK4 S3--------------!!
    bq%ur = bq%tOb(0)%p / bq%tOb(0)%tD
    bq%vr = bq%tOb(0)%q / bq%tOb(0)%tD
    bq%pbpr = bq%tOb(0)%p / bq%por
    bq%qbpr = bq%tOb(0)%q / bq%por   
    rTime=(bq%tOb(0)%rtm + bq%tOb(1)%rtm)/2d0
    call bq%dynaMatrices(rTime,bq%tOb(0)%tD, bq%ur, bq%vr)
    call bq%solveAll(rTime, bq%tOb(0)%p, bq%tOb(0)%q, &
      bq%pbpr, bq%qbpr, bq%presr, bq%tOb(0)%e, &
      bq%gXW, bq%gXE, bq%gXPQ, bq%gRE, bq%gRPQ, bq%sysC)
    call bq%updateSoln(3)
    !!-----------End RK4 S3------------!!

    !!-------------RK4 S4--------------!!
    bq%ur = bq%tOb(0)%p / bq%tOb(0)%tD
    bq%vr = bq%tOb(0)%q / bq%tOb(0)%tD
    bq%pbpr = bq%tOb(0)%p / bq%por
    bq%qbpr = bq%tOb(0)%q / bq%por   
    rTime=bq%tOb(0)%rtm
    call bq%dynaMatrices(rTime,bq%tOb(0)%tD, bq%ur, bq%vr)
    call bq%solveAll(rTime, bq%tOb(0)%p, bq%tOb(0)%q, &
      bq%pbpr, bq%qbpr, bq%presr, bq%tOb(0)%e, &
      bq%gXW, bq%gXE, bq%gXPQ, bq%gRE, bq%gRPQ, bq%sysC)
    call bq%updateSoln(4)
    !!-----------End RK4 S4------------!!
  
    call bq%postInstructs   


  enddo

  call system_clock(bq%sysC(2))
  close(bq%wpEle(-1))
  write(9,*)"[MSG] boussinesqQuad End"
  write(9,'(" [TIM] ",F15.4)')1d0*(bq%sysC(2)-bq%sysC(1))/bq%sysRate
  close(9)  
end program boussinesqQuad