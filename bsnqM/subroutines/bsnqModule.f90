module bsnqModule
use bsnqGlobVars
use waveFileModule
use outAbsModule
use shipMod
use meshFreeMod
use vertVelMod
implicit none

  interface
    subroutine paralution_init(nthreads) BIND(C)
      use, intrinsic :: ISO_C_BINDING, only : C_INT
      integer(kind=C_INT), value, intent(in)  :: nthreads
    end subroutine paralution_init

    subroutine paralution_stop() BIND(C)
    end subroutine paralution_stop        
  end interface
  
  private
  integer(kind=C_K1)::maxNePoi=30
  integer(kind=C_K1)::maxNeEle=10
  !! Number of elements in sysC and sysT
  integer(kind=C_K1),parameter::nSysC=10 

  type, public :: bsnqVars
    integer(kind=C_K1)::npl,npt
    real(kind=C_K2)::rtm
    real(kind=C_K2),allocatable::p(:),q(:),e(:),tD(:)
  contains
    procedure ::  initBsnqVars
  end type bsnqVars

  
  type, public :: bsnqCase
    
    character(len=C_KSTR)::probname,resumeFile
    integer(kind=C_K1)::npl,npq,npt,nele,nbnd,nbndtyp,nedg
    integer(kind=C_K1)::maxNePoi,maxNeEle,nbndp,nthrd,Sz(4)
    integer(kind=C_K1)::maxIter,fileOut,resumeOut    
    integer(kind=C_K1)::nnzl,nnzf,tStep,nTOb,nSOb
    integer(kind=C_K1),allocatable::conn(:,:),mabnd(:,:)    
    integer(kind=C_K1),allocatable::npoisur(:,:),bndP(:)
    integer(kind=C_K1),allocatable::bndPT(:)
    integer(kind=C_K1),allocatable::wpEle(:)
    integer(kind=C_K1),allocatable::ivl(:),linkl(:),ivq(:),linkq(:)            
    integer(kind=C_INT),allocatable::ivsl(:),jvsl(:)
    integer(kind=C_INT),allocatable::ivsf(:),jvsf(:)
    integer(kind=C_INT),allocatable::ele6x6(:,:),ele6x3(:,:)

    real(kind=C_K2)::dt,errLim,endTime,wvHReset
    real(kind=C_K2)::sysRate,sysT(nSysC)
    real(kind=C_K2),allocatable::cor(:,:),dep(:)    
    real(kind=C_K2),allocatable::invJ(:,:),bndS(:,:),bndPN(:,:)    
    real(kind=C_K2),allocatable::por(:),presr(:),vec6Tmp(:)
    real(kind=C_K2),allocatable::ur(:),vr(:),pbpr(:),qbpr(:)    
    real(kind=C_K2),allocatable::uhr(:),vhr(:)
    real(kind=C_K2),allocatable::mass1(:),mass2(:)    
    real(kind=C_K2),allocatable::rowMaxW(:),rowMaxE(:),rowMaxPQ(:)
    real(kind=C_K2),allocatable::gBs5(:),gBs6(:),absC(:)
    real(kind=C_K2),allocatable::gGx(:),gGy(:),gNAdv(:)    
    real(kind=C_K2),allocatable::gPGx(:),gPGy(:)
    real(kind=C_K2),allocatable::gCxF(:),gCyF(:),gDMat(:)
    real(kind=C_K2),allocatable::gBs1(:),gBs2(:),gBs3(:),gBs4(:)
    real(kind=C_K2),allocatable::etaMax(:),etaMin(:),wpLoc(:,:)
    real(kind=C_DOUBLE),allocatable::gXW(:),gXE(:),gXPQ(:)
    real(kind=C_DOUBLE),allocatable::gRE(:),gRPQ(:)
    real(kind=C_DOUBLE),allocatable::gMW(:),gME(:),gMPQ(:)

    !! Deallocated in destructR1
    integer(kind=C_K1),allocatable::p2p(:,:),p2e(:,:)
    integer(kind=C_K1),allocatable::ivf(:),linkf(:)  
    real(kind=C_K2),allocatable::massW(:),massE(:)
    real(kind=C_K2),allocatable::gFBs1(:),gFBs2(:),gFBs3(:),gFBs4(:)
    real(kind=C_K2),allocatable::aFull(:),gFW(:)

    integer(kind=C_KCLK)::sysC(nSysC)
    logical(kind=C_LG)::resume,presOn,absOn    
    type(wvFileType)::wvF
    type(shipType),allocatable::sh(:)
    type(absTyp),allocatable::absOb(:)
    type(bsnqVars),allocatable::tOb(:),sOb(:)

    !! Optional initialisation
    type(mfPoiTyp),allocatable::pObf(:)
    type(vertVelDerv),allocatable::bDf(:)
    

  contains    

    procedure ::  initMat
    procedure ::  meshRead
    procedure ::  femInit
    procedure ::  setRun
    procedure ::  statMatrices
    procedure ::  dynaMatrices
    procedure ::  CSRMatrices
    procedure ::  destructR1    
    procedure ::  outputXML    
    procedure ::  preInstructs
    procedure ::  postInstructs
    procedure ::  solveAll
    procedure ::  diriBCEtaDiff
    procedure ::  diriBCPQDiff
    procedure ::  diriBCEta
    procedure ::  diriBCPQ
    procedure ::  updateSoln
    procedure ::  writeWaveProbe
    procedure ::  getEtaPQForXY
    procedure ::  writeResume
    procedure ::  readResume
    !procedure ::  destructor
    
    procedure ::  setMFree
    procedure ::  calcVertVelDerv
    procedure ::  findEleForLocXY1    !For one location. No OpenMP
    procedure ::  findEleForLocXY2    !For a matrix of locs. OpenMP
    procedure ::  getVertVel          !using vertVelExp to calculate
    procedure ::  testGetVertVel          

  end type bsnqCase

contains

  include 'bsnqModuleFncs.f90'
  include 'bsnqModuleFncs2.f90'
  include 'initialCondition.f90'
  include 'outputXML.f90'

!!---------------------------preInstructs--------------------------!!
  subroutine preInstructs(b)
  implicit none

    class(bsnqCase),intent(inout)::b    
    integer(kind=C_K1)::i

    call system_clock(b%sysC(3)) 
    
    b%tStep=b%tStep+1            

    b%sysT(1)=0d0 ! To time PQ soln in Predictor + Corrector    
    b%sysT(2)=0d0 ! To time solveAll

    do i=b%nTOb-1,1,-1
      b%tOb(i)%rtm = b%tOb(i-1)%rtm
      b%tOb(i)%e = b%tOb(i-1)%e
      b%tOb(i)%p = b%tOb(i-1)%p
      b%tOb(i)%q = b%tOb(i-1)%q
      b%tOb(i)%tD = b%tOb(i-1)%tD
    enddo
    b%tOb(0)%rtm = b%tOb(1)%rtm + b%dt

    write(9,'(" ------Time : ",F20.6,"------")') b%tOb(0)%rtm

  end subroutine preInstructs
!!-------------------------End preInstructs------------------------!!



!!--------------------------postInstructs--------------------------!!
  subroutine postInstructs(b)
  implicit none

    class(bsnqCase),intent(inout)::b

    integer(kind=C_K1)::i    
    real(kind=C_K2)::tmpr1

    
    b%tOb(0)%e = b%tOb(1)%e + 1d0/6d0*(b%sOb(1)%e &
      + 2d0*b%sOb(2)%e + 2d0*b%sOb(3)%e + b%sOb(4)%e)
    b%tOb(0)%p = b%tOb(1)%p + 1d0/6d0*(b%sOb(1)%p &
      + 2d0*b%sOb(2)%p + 2d0*b%sOb(3)%p + b%sOb(4)%p)
    b%tOb(0)%q = b%tOb(1)%q + 1d0/6d0*(b%sOb(1)%q &
      + 2d0*b%sOb(2)%q + 2d0*b%sOb(3)%q + b%sOb(4)%q)

    !! Forcing Dirichlet BC
    call b%diriBCEta(b%tOb(0)%e,b%tOb(0)%rtm)
    call b%diriBCPQ(b%tOb(0)%p,b%tOb(0)%q,b%tOb(0)%rtm)
    
    b%tOb(0)%tD(1:b%npl) = b%dep(1:b%npl) + b%tOb(0)%e
    call fillMidPoiVals(b%npl,b%npt,b%nele,b%conn,b%tOb(0)%tD)

    do i=1,b%npl
      tmpr1=b%tOb(0)%e(i)
      if(b%etaMin(i).gt.tmpr1) b%etaMin(i)=tmpr1
      if(b%etaMax(i).lt.tmpr1) b%etaMax(i)=tmpr1
    enddo

    ! do i=1,b%nbndp
    !   j2=b%bndPT(i)
    !   if(j2.eq.11)then
    !     i2=b%bndP(i)
    !     call b%wvIn%getEta(b%tOb(0)%rtm,b%cor(i2,1),b%cor(i2,2),tmpr1)
    !     call b%wvIn%getEta(b%tOb(1)%rtm,b%cor(i2,1),b%cor(i2,2),tmpr2)
    !     write(9,302)"InWv",b%tOb(0)%rtm,b%tOb(0)%e(i2),&
    !       b%tOb(0)%p(i2),b%tOb(0)%q(i2)
    !     write(9,302)"InEt",b%tOb(1)%e(i2),tmpr1,tmpr2,tmpr1-tmpr2
    !     exit
    !   endif
    ! enddo

    if(allocated(b%bDf))then 
      call b%calcVertVelDerv    
      !call b%testGetVertVel
    endif

    if(mod(b%tStep,b%fileOut).eq.0) then
      call b%outputXML
    endif

    if(mod(b%tStep,b%resumeOut).eq.0) then
      call b%writeResume
    endif

    tmpr1=b%wvHReset
    if( mod( b%tOb(0)%rtm, tmpr1 ) .lt. 0.1*b%dt )then
      b%etaMin=b%tOb(0)%e
      b%etaMax=b%tOb(0)%e
    endif
    !! Ship drag calculation
    if(b%presOn)then
      if(b%sh(1)%dragFlag)then !Optional to calc drag
        b%vec6Tmp = b%tOb(0)%tD - b%dep
        i=b%npt
        call b%sh(1)%calcDrag(b%tOb(0)%rtm,i,b%cor(1:i,1),&
          b%cor(1:i,2),b%vec6Tmp(1:i),tmpr1)
        write(9,303)'shFx',b%tOb(0)%rtm,tmpr1
      endif
    endif

    !write(201,*)

    call b%writeWaveProbe    

    call system_clock(b%sysC(4))
    ! bq%sysT(1) = To time PQ soln in Predictor + Corrector    
    ! bq%sysT(2) = To time solveAll    
    ! bq%sysT(3) = To time the current time-step
    b%sysT(3)=1d0*(b%sysC(4)-b%sysC(3))/b%sysRate    
    write(9,301)"[TIL]", b%sysT(3), b%sysT(1)/b%sysT(3)*100d0, &
      b%sysT(2)/b%sysT(3)*100d0
    write(9,*)
    301 format('      |',a6,3F15.4)
    302 format('      |',a6,4F15.4)
    303 format('      |',a6,2F15.4)

  end subroutine postInstructs
!!------------------------End postInstructs------------------------!!



!!----------------------------updateSoln---------------------------!!
  subroutine updateSoln(b,step)
  implicit none

    class(bsnqCase),intent(inout)::b    
    integer(kind=4),intent(in)::step

    b%sysT(1)=b%sysT(1)+1d0*(b%sysC(8)-b%sysC(7))/b%sysRate
    b%sysT(2)=b%sysT(2)+1d0*(b%sysC(6)-b%sysC(5))/b%sysRate
    
    select case(step)

      case (1)
        b%sOb(1)%e = b%gXE
        b%sOb(1)%p = b%gXPQ(1:b%npt)
        b%sOb(1)%q = b%gXPQ(b%npt+1:2*b%npt)        
        
        b%tOb(0)%e = b%tOb(1)%e + b%gXE/2d0
        b%tOb(0)%p = b%tOb(1)%p + b%gXPQ(1:b%npt)/2d0
        b%tOb(0)%q = b%tOb(1)%q + b%gXPQ(b%npt+1:2*b%npt)/2d0
        b%tOb(0)%tD(1:b%npl) = b%dep(1:b%npl) + b%tOb(0)%e
        call fillMidPoiVals(b%npl,b%npt,b%nele,b%conn,b%tOb(0)%tD) 

      case (2)
        b%sOb(2)%e = b%gXE
        b%sOb(2)%p = b%gXPQ(1:b%npt)
        b%sOb(2)%q = b%gXPQ(b%npt+1:2*b%npt)        
        
        b%tOb(0)%e = b%tOb(1)%e + b%gXE/2d0
        b%tOb(0)%p = b%tOb(1)%p + b%gXPQ(1:b%npt)/2d0
        b%tOb(0)%q = b%tOb(1)%q + b%gXPQ(b%npt+1:2*b%npt)/2d0
        b%tOb(0)%tD(1:b%npl) = b%dep(1:b%npl) + b%tOb(0)%e
        call fillMidPoiVals(b%npl,b%npt,b%nele,b%conn,b%tOb(0)%tD) 

      case (3)
        b%sOb(3)%e = b%gXE
        b%sOb(3)%p = b%gXPQ(1:b%npt)
        b%sOb(3)%q = b%gXPQ(b%npt+1:2*b%npt)        
        
        b%tOb(0)%e = b%tOb(1)%e + b%gXE
        b%tOb(0)%p = b%tOb(1)%p + b%gXPQ(1:b%npt)
        b%tOb(0)%q = b%tOb(1)%q + b%gXPQ(b%npt+1:2*b%npt)
        b%tOb(0)%tD(1:b%npl) = b%dep(1:b%npl) + b%tOb(0)%e
        call fillMidPoiVals(b%npl,b%npt,b%nele,b%conn,b%tOb(0)%tD) 

      case (4)
        b%sOb(4)%e = b%gXE
        b%sOb(4)%p = b%gXPQ(1:b%npt)
        b%sOb(4)%q = b%gXPQ(b%npt+1:2*b%npt)        

    end select        

  end subroutine updateSoln
!!--------------------------End updateSoln-------------------------!!



!!----------------------------diriBCEta----------------------------!!
  subroutine diriBCEta(b,mat,rTt0)
  implicit none

    class(bsnqCase),intent(in)::b    
    real(kind=C_K2),intent(in)::rTt0
    real(kind=C_K2),intent(inout)::mat(b%npl)    
    integer(kind=C_K1)::i,i2,j2
    real(kind=C_K2)::tmpr1

    !! Note : Consistent with SemiDirect only
    call b%wvF%getEta(rTt0,tmpr1)    
    do i=1,b%nbndp
      i2=b%bndP(i)
      j2=b%bndPT(i)
      if(i2.gt.b%npl)cycle !Linear only
      if(j2.eq.11)then        
        mat(i2)=tmpr1

      elseif(j2.eq.14)then
        mat(i2)=0d0

      endif
    enddo
  end subroutine diriBCEta
!!--------------------------End diriBCEta--------------------------!!



!!----------------------------diriBCPQ-----------------------------!!
  subroutine diriBCPQ(b,p,q,rTt0)
  implicit none

    class(bsnqCase),intent(in)::b    
    real(kind=C_K2),intent(in)::rTt0
    real(kind=C_K2),intent(inout)::p(b%npt),q(b%npt)
    integer(kind=C_K1)::i,i2,j2
    real(kind=C_K2)::tmpr1,tmpr2

    !! Note : Consistent with SemiDirect only
    call b%wvF%getPQ(rTt0,tmpr1,tmpr2)    
    do i=1,b%nbndp
      i2=b%bndP(i)
      j2=b%bndPT(i)
      if(j2.eq.11)then        
        p(i2)=tmpr1
        q(i2)=tmpr2

      elseif((j2.eq.12).or.(j2.eq.14))then
        p(i2)=0d0
        q(i2)=0d0

      elseif(j2.eq.13)then
        if(abs(b%bndPN(i2,1)).gt.0.1)then
          p(i2)=0d0
        endif
        if(abs(b%bndPN(i2,2)).gt.0.1)then
          q(i2)=0d0
        endif

      endif
    enddo
  end subroutine diriBCPQ
!!--------------------------End diriBCPQ---------------------------!!



!!--------------------------diriBCEtaDiff--------------------------!!
  subroutine diriBCEtaDiff(b,mat,rTt0,rTt1)
  implicit none

    class(bsnqCase),intent(in)::b
    real(kind=C_K2),intent(in)::rTt0,rTt1
    real(kind=C_K2),intent(inout)::mat(b%npl)    
    integer(kind=C_K1)::i,i2,j2
    real(kind=C_K2)::tmpr1,tmpr2

    !! Note : Consistent with SemiDirect only
    call b%wvF%getEta(rTt0,tmpr1)
    call b%wvF%getEta(rTt1,tmpr2)
    !write(201,'(3F15.6)',advance='no')rTt0,tmpr1,tmpr2
    do i=1,b%nbndp
      i2=b%bndP(i)
      j2=b%bndPT(i)
      if(i2.gt.b%npl)cycle !Linear only
      if(j2.eq.11)then        
        mat(i2)=tmpr1-tmpr2

      elseif(j2.eq.14)then
        mat(i2)=0d0

      endif
    enddo
  end subroutine diriBCEtaDiff
!!------------------------End diriBCEtaDiff------------------------!!



!!--------------------------diriBCPQDiff---------------------------!!
  subroutine diriBCPQDiff(b,mat,rTt0,rTt1)
  implicit none

    class(bsnqCase),intent(in)::b
    real(kind=C_K2),intent(in)::rTt0,rTt1
    real(kind=C_K2),intent(inout)::mat(2*b%npt)
    integer(kind=C_K1)::i,i2,j2
    real(kind=C_K2)::tmpr1,tmpr2,tmpr3,tmpr4

    !! Note : Consistent with SemiDirect only
    call b%wvF%getPQ(rTt0,tmpr1,tmpr2)         
    call b%wvF%getPQ(rTt1,tmpr3,tmpr4)
    !write(201,'(2F15.6)',advance='no')tmpr1,tmpr3
    do i=1,b%nbndp
      i2=b%bndP(i)
      j2=b%bndPT(i)
      if(j2.eq.11)then        
        mat(i2)=tmpr1-tmpr3
        mat(b%npt+i2)=tmpr2-tmpr4

      elseif((j2.eq.12).or.(j2.eq.14))then
        mat(i2)=0d0
        mat(b%npt+i2)=0d0

      elseif(j2.eq.13)then
        if(abs(b%bndPN(i2,1)).gt.0.1)then
          mat(i2)=0d0
        endif
        if(abs(b%bndPN(i2,2)).gt.0.1)then
          mat(b%npt+i2)=0d0
        endif

      endif
    enddo
  end subroutine diriBCPQDiff
!!------------------------End diriBCPQDiff-------------------------!!



!!-----------------------------initMat-----------------------------!!
  subroutine initMat(b)
  implicit none

    class(bsnqCase),intent(inout)::b
    integer(kind=C_K1)::i,i1,j,j1

    i=b%npl
    j=b%npt
    i1=b%ivl(0)
    j1=b%ivq(0)

    !Allocations
    allocate(b%gXE(i),b%gXPQ(2*j),b%gXW(i))
    allocate(b%por(j),b%vec6Tmp(j))
    allocate(b%ur(j),b%vr(j),b%pbpr(j),b%qbpr(j))
    allocate(b%uhr(j),b%vhr(j))

    allocate(b%massW(i1*i),b%massE(i1*i))
    allocate(b%mass1(j1*j),b%mass2(i1*i))    
    allocate(b%gBs1(j1*j),b%gBs2(j1*j))
    allocate(b%gBs3(j1*j),b%gBs4(j1*j))
    allocate(b%gCxF(j1*i),b%gCyF(j1*i),b%gDMat(i1*i))
    allocate(b%gBs5(i1*j),b%gBs6(i1*j),b%absC(j))
    allocate(b%gGx(i1*j),b%gGy(i1*j),b%gNAdv(j1*j))
    allocate(b%gPGx(j1*j),b%gPGy(j1*j))
    allocate(b%gFBs1(j1*j),b%gFBs2(j1*j))
    allocate(b%gFBs3(j1*j),b%gFBs4(j1*j),b%gFW(i1*i))
    allocate(b%aFull(b%ivf(0)*2*j))
    allocate(b%rowMaxW(i),b%rowMaxE(i),b%rowMaxPQ(2*j))
    allocate(b%gRE(i),b%gRPQ(2*j))
    allocate(b%etaMax(i),b%etaMin(i),b%presr(j))
    allocate(b%ele6x6(b%nele,36),b%ele6x3(b%nele,18))

    b%Sz(1)=i1*i ![3x3] ![ivl(0) * npl]
    b%Sz(2)=j1*i ![3x6] ![ivq(0) * npl]
    b%Sz(3)=i1*j ![6x3] ![ivl(0) * npt]
    b%Sz(4)=j1*j ![6x6] ![ivq(0) * npt]
    b%por=1d0
    
    b%absC=0d0
    if(b%absOn)then
      do i=1,b%absOb(1)%N
        call b%absOb(i)%calcAbsC(b%npt,b%cor,b%absC)
      enddo  
    endif

    b%nTOb=2
    allocate(b%tOb(0:b%nTOb-1))
    do i=0,b%nTOb-1
      call b%tOb(i)%initBsnqVars(b%npl,b%npt)
    enddo

    b%nSOb=4
    allocate(b%sOb(b%nSOb))
    do i=1,b%nSOb
      call b%sOb(i)%initBsnqVars(b%npl,b%npt)
    enddo    

    b%tStep=0

    if(b%resume)then
      ! [Note]:
      ! Check bsnqM/Dev - Logs/log_bsnqM_v0003.md
      ! to understand why the variables 
      ! gXW, gXE, gXPQ are being written in 
      ! resume file.
      call b%readResume
    
    else
      b%tOb(0)%rtm=0d0
      b%tOb(0)%e=0d0
      b%tOb(0)%p=0d0
      b%tOb(0)%q=0d0      

      b%gXW = 0d0
      b%gXE = 0d0
      b%gXPQ = 0d0
    endif
    
    b%tOb(0)%tD(1:b%npl)=b%dep(1:b%npl)+b%tOb(0)%e
    call fillMidPoiVals(b%npl,b%npt,b%nele,b%conn,b%tOb(0)%tD)    

    ! call solitIC(b%npl, b%npt, b%cor, &
    !   b%tOb(0)%e, b%tOb(0)%p, b%tOb(0)%q)
    
    b%etaMin=b%tOb(0)%e
    b%etaMax=b%tOb(0)%e    

    b%presr=0d0    

    ! call testMls2DDx
    ! stop

    call paralution_init(b%nthrd)    

    call b%outputXML

    write(9,*)"[MSG] Done initMat"
    write(9,*)

  end subroutine initMat
!!---------------------------End initMat---------------------------!!



!!-------------------------writeWaveProbe--------------------------!!
  subroutine writeWaveProbe(b)
  implicit none

    class(bsnqCase),intent(in)::b
    integer(kind=C_K1)::nq(6),i,i2,k
    real(kind=C_K2)::N3(3),N6(6),tmpr1,tmpr2,tmpr3

    !Writing probes value
    k=b%wpEle(-1)
    write(k,'(F15.6)',advance='no')b%tOb(0)%rtm
    do i=1,b%wpEle(0)

      tmpr1=0d0
      tmpr2=0d0
      tmpr3=0d0

      if(b%wpEle(i).ne.-1)then

        nq=b%conn(b%wpEle(i),:)
        call fem_N6i(b%wpLoc(i,3), b%wpLoc(i,4), N6)
        call fem_N3i(b%wpLoc(i,3), b%wpLoc(i,4), N3)

        do i2=1,3
          tmpr1 = tmpr1 + b%tOb(0)%e(nq(i2)) * N3(i2) 
        enddo
        do i2=1,6
          tmpr2 = tmpr2 + b%tOb(0)%p(nq(i2)) * N6(i2) 
          tmpr3 = tmpr3 + b%tOb(0)%q(nq(i2)) * N6(i2)                     
        enddo

      endif

      write(k,'(5F15.6)',advance='no')b%wpLoc(i,1), b%wpLoc(i,2),&
        tmpr1, tmpr2, tmpr3
    enddo
    write(k,*)

  end subroutine writeWaveProbe
!!-----------------------End writeWaveProbe------------------------!!



!!--------------------------getEtaPQForXY--------------------------!!
  subroutine getEtaPQForXY(b,np,xin,yin,eta,p,q,wrki,wrkr,err)
  implicit none

    class(bsnqCase),intent(in)::b    
    integer(kind=C_K1),intent(in)::np
    real(kind=C_K2),intent(in)::xin(np),yin(np)
    integer(kind=C_K1),intent(out)::wrki(np),err(np)
    real(kind=C_K2),intent(out)::eta(np),p(np),q(np),wrkr(np,2)

    integer(kind=C_K1)::nq(6),i,k
    real(kind=C_K2)::wei(6),etaLoc,pLoc,qLoc,hLoc

    call b%findEleForLocXY2(np,xin,yin,wrki,wrkr(:,1),wrkr(:,2))

    !$OMP PARALLEL DEFAULT(shared) &
    !$OMP   PRIVATE(i, nq, wei, hLoc, etaLoc, pLoc, qLoc, k)
    !$OMP DO SCHEDULE(dynamic,10)
    do i=1,np

      if(wrki(i).eq.-1)then
        err(i)=1
        eta(i)=0d0
        p(i)=0d0
        q(i)=0d0   
        cycle     
      endif

      nq = b%conn(wrki(i),:)
      call fem_N6i(wrkr(i,1),wrkr(i,2),wei)

      hLoc=0d0
      etaLoc=0d0
      pLoc=0d0
      qLoc=0d0

      do k=1,6
        hLoc = hLoc + wei(k) * b%dep(nq(k))
        etaLoc = etaLoc + wei(k) * b%tOb(0)%tD(nq(k)) !Using totalDep 
        pLoc = pLoc + wei(k) * b%tOb(0)%p(nq(k))
        qLoc = qLoc + wei(k) * b%tOb(0)%q(nq(k))
      enddo
      etaLoc = etaLoc - hLoc

      err(i)=0
      eta(i)=etaLoc
      p(i)=pLoc
      q(i)=qLoc

    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL    

  end subroutine getEtaPQForXY

!!------------------------End getEtaPQForXY------------------------!!



!!--------------------------initBsnqVars---------------------------!!
  subroutine initBsnqVars(b,npl,npt)
  use bsnqGlobVars
  implicit none

    class(bsnqVars),intent(inout)::b
    integer(kind=C_K1),intent(in)::npl,npt
    logical(kind=C_LG)::tmp

    tmp=allocated(b%e).or.allocated(b%p).or.&
      allocated(b%q).or.allocated(b%tD)    
    if(tmp)then
      write(9,'(" [ERR] Already allocated bsnqVars")')
      stop
    endif
    allocate(b%e(npl),b%p(npt),b%q(npt),b%tD(npt))
    b%npl=npl
    b%npt=npt
    b%e=0d0
    b%p=0d0
    b%q=0d0
    b%tD=0d0

  end subroutine initBsnqVars
!!------------------------End initBsnqVars-------------------------!!


end module bsnqModule