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
!! Time-Stepping : AdBs3E : Adam-Bashforth Explicit 3-point

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

    !!------------AdBs3E S1------------!!
    bq%ur = bq%tOb(1)%p / bq%tOb(1)%tD
    bq%vr = bq%tOb(1)%q / bq%tOb(1)%tD
    bq%pbpr = bq%tOb(1)%p / bq%por
    bq%qbpr = bq%tOb(1)%q / bq%por   
    rTime=bq%tOb(1)%rtm
    call bq%dynaMatrices(bq%tOb(1)%tD, bq%ur, bq%vr)
    call bq%solveAll(rTime, bq%tOb(1)%p, bq%tOb(1)%q, &
      bq%pbpr, bq%qbpr, bq%tOb(1)%e, &
      bq%gXW, bq%gXE, bq%gXPQ, bq%gRE, bq%gRPQ)
    call bq%updateSoln
    !!----------End AdBs3E S1----------!!

    call bq%postInstructs

  enddo

  call system_clock(bq%sysC(2))
  write(9,*)"[MSG] boussinesqQuad End"
  write(9,'(" [TIM] ",F15.4)')1d0*(bq%sysC(2)-bq%sysC(1))/bq%sysRate
  close(9)
end program boussinesqQuad