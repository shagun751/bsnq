module bsnqModule
use bsnqGlobVars
use waveFileModule
use outAbsModule
use shipMod
use meshFreeMod
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
  integer(kind=C_K1)::i,i1,i2,j,j1,j2,k,k1,k2,mf
  integer(kind=C_K1)::iel,n1,n2,n3,n4,n5,n6
  integer(kind=C_K1)::nq(6),nl(3)
  integer(kind=C_K1)::tmpi1,tmpi2,tmpi3,tmpi4,tmpi5
  integer(kind=C_K1)::maxNePoi=30
  integer(kind=C_K1)::maxNeEle=10

  integer(kind=8)::sysC(10)
  real(kind=C_K2)::tmpr1,tmpr2,tmpr3,tmpr4,tmpr5,tmpr6
  character(len=C_KSTR)::bqtxt  
  logical(kind=C_LG)::ex  

  type, public :: bsnqVars
    integer(kind=C_K1)::npl,npt
    real(kind=C_K2)::rtm
    real(kind=C_K2),allocatable::p(:),q(:),e(:),tD(:)
  contains
    procedure ::  initBsnqVars
  end type bsnqVars

  type, public :: bsnqDerv
    integer(kind=C_K1)::np
    real(kind=C_K2)::rtm
    real(kind=C_K2),allocatable::ux(:),uxx(:),uxxx(:)
    real(kind=C_K2),allocatable::px(:),pxx(:),pxxx(:)
  contains
    procedure ::  init => initBsnqDerv
  end type bsnqDerv

  
  type, public :: bsnqCase
    
    character(len=C_KSTR)::probname,resumeFile
    integer(kind=C_K1)::npl,npq,npt,nele,nbnd,nbndtyp,nedg
    integer(kind=C_K1)::maxNePoi,maxNeEle,nbndp,nthrd,Sz(4)
    integer(kind=C_K1)::maxIter,fileOut,resumeOut    
    integer(kind=C_K1)::nnzl,nnzf,tStep,nTOb,nSOb
    integer(kind=C_K1),allocatable::conn(:,:),mabnd(:,:)    
    integer(kind=C_K1),allocatable::npoisur(:,:),bndP(:)
    integer(kind=C_K1),allocatable::bndPT(:),probe(:)
    integer(kind=C_K1),allocatable::ivl(:),linkl(:),ivq(:),linkq(:)            
    integer(kind=C_INT),allocatable::ivsl(:),jvsl(:)
    integer(kind=C_INT),allocatable::ivsf(:),jvsf(:)
    integer(kind=C_INT),allocatable::ele6x6(:,:),ele6x3(:,:)

    real(kind=C_K2)::dt,errLim,endTime,rTime,wvHReset
    real(kind=C_K2)::sysRate,sysT(10)
    real(kind=C_K2),allocatable::cor(:,:),dep(:)    
    real(kind=C_K2),allocatable::invJ(:,:),bndS(:,:),bndPN(:,:)
    real(kind=C_K2),allocatable::por(:),presr(:),vec6Tmp(:)
    real(kind=C_K2),allocatable::ur(:),vr(:),pbpr(:),qbpr(:)    
    real(kind=C_K2),allocatable::mass1(:),mass2(:)    
    real(kind=C_K2),allocatable::rowMaxW(:),rowMaxE(:),rowMaxPQ(:)
    real(kind=C_K2),allocatable::gBs5(:),gBs6(:),absC(:)
    real(kind=C_K2),allocatable::gGx(:),gGy(:),gNAdv(:)    
    real(kind=C_K2),allocatable::gPGx(:),gPGy(:)
    real(kind=C_K2),allocatable::gCxF(:),gCyF(:),gDMat(:)
    real(kind=C_K2),allocatable::gBs1(:),gBs2(:),gBs3(:),gBs4(:)
    real(kind=C_K2),allocatable::etaMax(:),etaMin(:),probeLoc(:,:)
    real(kind=C_DOUBLE),allocatable::gXW(:),gXE(:),gXPQ(:)
    real(kind=C_DOUBLE),allocatable::gRE(:),gRPQ(:)
    real(kind=C_DOUBLE),allocatable::gMW(:),gME(:),gMPQ(:)

    !! Deallocated in destructR1
    integer(kind=C_K1),allocatable::p2p(:,:),p2e(:,:)
    integer(kind=C_K1),allocatable::ivf(:),linkf(:)  
    real(kind=C_K2),allocatable::massW(:),massE(:)
    real(kind=C_K2),allocatable::gFBs1(:),gFBs2(:),gFBs3(:),gFBs4(:)
    real(kind=C_K2),allocatable::aFull(:),gFW(:)

    integer(kind=8)::sysC(10)
    logical(kind=C_LG)::resume,presOn,absOn    
    type(wvFileType)::wvF
    type(shipType),allocatable::sh(:)
    type(absTyp),allocatable::absOb(:)
    type(bsnqVars),allocatable::tOb(:),sOb(:)

    !! Optional initialisation
    type(mfPoiTyp),allocatable::pObf(:)
    type(bsnqDerv)::bDf
    

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
    !procedure ::  destructor
    procedure ::  setMFree
    procedure ::  calcDerv

  end type bsnqCase

contains

  include 'bsnqModuleFncs.f90'
  include 'bsnqModuleFncs2.f90'
  include 'outputXML.f90'

!!---------------------------preInstructs--------------------------!!
  subroutine preInstructs(b)
  implicit none

    class(bsnqCase),intent(inout)::b    

    b%tStep=b%tStep+1    
    b%rTime=b%rTime+b%dt
    write(9,'(" ------Time : ",I10," ",F20.6,"------")')&
      b%tStep,b%rTime
    call system_clock(b%sysC(3)) 

    b%sysT(1)=0d0 ! To time PQ soln in Predictor + Corrector    
    b%sysT(2)=0d0 ! solveAll

    do i=b%nTOb-1,1,-1
      b%tOb(i)%rtm = b%tOb(i-1)%rtm
      b%tOb(i)%e = b%tOb(i-1)%e
      b%tOb(i)%p = b%tOb(i-1)%p
      b%tOb(i)%q = b%tOb(i-1)%q
      b%tOb(i)%tD = b%tOb(i-1)%tD
    enddo
    b%tOb(0)%rtm=b%rTime

  end subroutine preInstructs
!!-------------------------End preInstructs------------------------!!



!!--------------------------postInstructs--------------------------!!
  subroutine postInstructs(b)
  implicit none

    class(bsnqCase),intent(inout)::b    

    
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

    if(allocated(b%pObf)) call b%calcDerv

    if(mod(b%tStep,b%fileOut).eq.0) then
      call b%outputXML
    endif

    tmpr1=b%wvHReset
    if(mod(b%rTime,tmpr1).lt.0.05*tmpr1)then
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

    !Writing probes values
    k=b%probe(-1)
    write(k,'(F15.6)',advance='no')b%tOb(0)%rtm
    do i=1,b%probe(0)
      j=b%probe(i)
      write(k,'(5F15.6)',advance='no')b%cor(j,1), b%cor(j,2),&
        b%tOb(0)%e(j), b%tOb(0)%p(j), b%tOb(0)%q(j)
    enddo
    write(k,*)

    call system_clock(b%sysC(4))
    tmpr1=1d0*(b%sysC(4)-b%sysC(3))/b%sysRate
    tmpr2=b%sysT(1)
    !write(9,301)"[TIL]",tmpr1,tmpr2,tmpr2/tmpr1*100d0
    write(9,301)"[TIL]",tmpr1,tmpr2/tmpr1*100d0,&
      b%sysT(2)/tmpr1*100d0
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

    b%sysT(1)=b%sysT(1)+1d0*(sysC(2)-sysC(1))/b%sysRate
    b%sysT(2)=b%sysT(2)+1d0*(sysC(4)-sysC(3))/b%sysRate
    
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



!!-----------------------------solitIC-----------------------------!!
  subroutine solitIC(npl,npt,cor,er,pr,qr)
  use bsnqGlobVars
  implicit none
    
    integer(kind=C_K1),intent(in)::npl,npt
    real(kind=C_K2),intent(in)::cor(npt,2)
    real(kind=C_K2),intent(out)::er(npl),pr(npt),qr(npt)

    ! pr=0d0
    ! qr=0d0
    ! !er=0.045d0*dexp(-2d0*( (cor(1:npl,1)-18.288d0)**2 ))
    ! do i=1,npl
    !   er(i)=(cor(i,1)-18.288d0)**2
    !   if(er(i).gt.25)then
    !     er(i)=0d0
    !   else
    !     er(i)=0.045*exp(-2d0*er(i))
    !   endif
    ! enddo

    er=0d0
    pr=0d0
    qr=0d0
    tmpr1=0.45d0
    tmpr2=0.045d0
    tmpr3=dsqrt(grav*(tmpr1+tmpr2))
    tmpr4=dsqrt(3*tmpr2/(4*(tmpr1**3)))
    do i=1,npt
      if((cor(i,1).ge.3d0).and.(cor(i,1).le.19d0)) then
        tmpr5=tmpr2/(dcosh(tmpr4*(cor(i,1)-(11d0)))**2)
        pr(i)=tmpr3*tmpr5
        if(i.le.npl) then
          er(i)=tmpr5
        endif
      endif
    enddo  

  end subroutine solitIC
!!---------------------------End solitIC---------------------------!!



!!----------------------------diriBCEta----------------------------!!
  subroutine diriBCEta(b,mat,rTt0)
  implicit none

    class(bsnqCase),intent(in)::b    
    real(kind=C_K2),intent(in)::rTt0
    real(kind=C_K2),intent(inout)::mat(b%npl)    

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



!!--------------------------initBsnqDerv---------------------------!!
  subroutine initBsnqDerv(b,np)
  use bsnqGlobVars
  implicit none

    class(bsnqDerv),intent(inout)::b
    integer(kind=C_K1),intent(in)::np

    allocate(b%ux(np), b%uxx(np), b%uxxx(np))
    allocate(b%px(np), b%pxx(np), b%pxxx(np))
    b%np=np    

  end subroutine initBsnqDerv
!!------------------------End initBsnqDerv-------------------------!!

end module bsnqModule