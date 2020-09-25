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
  type(bsnqCase)::bq  
  character(len=C_KSTR)::bqtxt

  real(kind=C_K2),allocatable::u(:),v(:),eta(:)
  real(kind=C_K2),allocatable::sedNelX(:),sedNelY(:)
  real(kind=C_K2),allocatable::sedvanX(:),sedVanY(:)
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

  allocate( u(bq%npt), v(bq%npt), eta(bq%npt) )
  allocate( sedNelX(bq%npt), sedNelY(bq%npt) )
  allocate( sedVanX(bq%npt), sedVanY(bq%npt) )

  call bq%caseOutputs
  
  do while(abs(bq%tOb(0)%rtm-bq%endTime).gt.bq%dt/2d0)
  
    call bq%preInstructs

    call bq%timeStepRK4
  
    call bq%postInstructs   

    u = bq%tOb(0)%p / bq%tOb(0)%tD
    v = bq%tOb(0)%q / bq%tOb(0)%tD
    eta = bq%tOb(0)%tD - bq%dep 

    call sedimentTrans(bq%npt, bq%dep, bq%cor(:,1), bq%cor(:,2), &
      u, v, eta, sedNelX, sedNelY, sedVanX, sedVanY)

    bq%sedNelX = sedNelX
    bq%sedNelY = sedNelY
    bq%sedVanX = sedVanX
    bq%sedVanY = sedVanY

    call bq%caseOutputs

  enddo

  call system_clock(bq%sysC(2))
  close(bq%wpEle(-1))
  write(9,*)"[MSG] boussinesqQuad End"
  write(9,'(" [TIM] ",F15.4)')1d0*(bq%sysC(2)-bq%sysC(1))/bq%sysRate
  close(9)  
end program boussinesqQuad



subroutine sedimentTrans(npt, dep, corx, cory, u, v, eta, &
  sedNelX, sedNelY, sedVanX, sedVanY)
use bsnqGlobVars
implicit none

  integer(kind=C_K1),intent(in):: npt
  real(kind=C_K2),intent(in):: dep(npt), u(npt), v(npt), eta(npt)
  real(kind=C_K2),intent(in):: corx(npt), cory(npt)
  real(kind=C_K2),intent(out):: sedNelX(npt), sedNelY(npt)
  real(kind=C_K2),intent(out):: sedVanX(npt), sedVanY(npt)

  !! Note: 
  !! Do not allocate the above variables again.
  !! They have already been allocated

  integer(kind=C_K1):: i
  real(kind=C_K2):: tmpr  

  sedNelX = 0d0
  sedNelY = 0d0
  sedVanX = 0d0
  sedVanY = 0d0

end subroutine sedimentTrans