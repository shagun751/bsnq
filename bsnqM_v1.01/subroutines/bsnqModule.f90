module bsnqModule
use bsnqGlobVars
implicit none
  
  private
  integer(kind=C_K1)::i,i1,i2,j,j1,j2,k,k1,k2
  integer(kind=C_K1)::iel,n1,n2,n3,n4,n5,n6
  integer(kind=C_K1)::nq(6),nl(3)
  integer(kind=C_K1)::tmpi1,tmpi2,tmpi3,tmpi4,tmpi5
  integer(kind=C_K1)::maxNePoi=30
  integer(kind=C_K1)::maxNeEle=10

  real(kind=C_K2)::tmpr1,tmpr2,tmpr3,tmpr4,tmpr5 
  character(len=C_KSTR)::bqtxt  
  logical::ex

  
  type, public :: waveType    
    real(kind=C_K2)::T,d,H,L,k,w
    real(kind=C_K2)::x0,y0
    real(kind=C_K2)::thDeg,thRad,csth,snth      
  end type waveType

  interface waveType
     procedure :: waveLenCalc
  end interface waveType

  
  type, public :: bsnqCase
    character(len=C_KSTR)::probname,resumeFile
    integer(kind=C_K1)::npl,npq,npt,nele,nbnd,nbndtyp,nedg
    integer(kind=C_K1)::maxNePoi,maxNeEle,nbndp,nthrd
    integer(kind=C_K1)::nTSteps,maxIter,fileOut,resumeOut    
    integer(kind=C_K1),allocatable::conn(:,:),mabnd(:,:)
    integer(kind=C_K1),allocatable::p2p(:,:),p2e(:,:)
    integer(kind=C_K1),allocatable::npoisur(:,:),bndP(:)
    integer(kind=C_K1),allocatable::bndPT(:)
    integer(kind=C_K1),allocatable::ivl(:),linkl(:),ivq(:),linkq(:)        
    integer(kind=C_K1),allocatable::ivf(:),linkf(:)    

    real(kind=C_K2)::dt,errLim,rTime
    real(kind=C_K2)::sysRate,sysTime(10)
    real(kind=C_K2),allocatable::cor(:,:),dep(:)    
    real(kind=C_K2),allocatable::invJ(:,:),bndS(:,:),bndPN(:,:)
    real(kind=C_K2),allocatable::pt0(:),qt0(:),et0(:),tDt0(:)
    real(kind=C_K2),allocatable::pt1(:),qt1(:),et1(:),tDt1(:)
    real(kind=C_K2),allocatable::pr(:),qr(:),er(:),tDr(:)
    real(kind=C_K2),allocatable::wr(:),por(:)
    real(kind=C_K2),allocatable::ur(:),vr(:),pbpr(:),qbpr(:)
    real(kind=C_K2),allocatable::massW(:),massE(:),massP(:),massQ(:)
    real(kind=C_K2),allocatable::gBs1(:),gBs2(:),gBs3(:),gBs4(:)
    real(kind=C_K2),allocatable::gCx(:),gCy(:),gDMat(:)
    real(kind=C_K2),allocatable::gBs5(:),gBs6(:)
    real(kind=C_K2),allocatable::gGx(:),gGy(:),gNAdv(:)
    real(kind=C_K2),allocatable::gFBs1(:),gFBs2(:),gFBs3(:),gFBs4(:)


    logical::resume,presOn

    integer(kind=8)::sysClk(10)
    integer(kind=C_INT),allocatable::ivCSR1(:),jvCSR1(:)
    integer(kind=C_INT),allocatable::ivCSR2(:),jvCSR2(:)

  contains    

    procedure ::  initMat
    procedure ::  meshRead
    procedure ::  femInit
    procedure ::  setRun
    procedure ::  statMatrices
    procedure ::  dynaMatrices
    procedure ::  destructR1
    !procedure ::  destructor

  end type bsnqCase

contains

!!---------------------------waveLenCalc---------------------------!!
type(waveType) function waveLenCalc(inT,inD,inH,inX0,inY0,inThDeg)
implicit none

  !! WaveLen using dispersion relation from Airy wave-theory

  integer(kind=C_K1)::iterMax
  real(kind=C_K2)::l0,newl,oldl,x,errLim
  real(kind=C_K2),intent(in)::inT,inD,inH,inX0,inY0,inThDeg


  waveLenCalc%T=inT
  waveLenCalc%d=inD
  waveLenCalc%H=inH
  waveLenCalc%x0=inX0
  waveLenCalc%y0=inY0
  waveLenCalc%thDeg=inThDeg
  waveLenCalc%w=2d0*pi/waveLenCalc%T
  waveLenCalc%thRad=waveLenCalc%thDeg*pi/180d0
  waveLenCalc%csth=dcos(waveLenCalc%thRad)
  waveLenCalc%snth=dsin(waveLenCalc%thRad)

  iterMax=50000
  errLim=1d-6
  l0 = (grav/2d0/pi)*(waveLenCalc%T)**2
  oldl = l0  
  do i = 1,iterMax
    newl = l0*dtanh(2d0*pi*(waveLenCalc%d)/oldl)
    !if(mod(i,100).eq.0) write(*,*)i,oldl,newl    
    x = abs(newl-oldl)
    if (x.le.errLim) then
      waveLenCalc%L = newl
      exit
    else
      oldl = newl
    end if
  end do

  if(i.ge.iterMax) then
    write(*,*)
    write(*,*)"[ERR] waveCalculator Error waveL",waveLenCalc%L
    write(*,*)
    waveLenCalc%L=-999
    waveLenCalc%k=0d0
    return
  endif

  waveLenCalc%k=2d0*pi/waveLenCalc%L

end function waveLenCalc
!!-------------------------End waveLenCalc-------------------------!!



!!-----------------------------meshRead----------------------------!!
  subroutine meshRead(b)
  implicit none

    class(bsnqCase),intent(inout)::b
    integer(kind=C_k1)::mf

    !Opening mesh file  
    bqtxt=trim(b%probname)//'.plt'
    inquire(file=trim(bqtxt),exist=ex)
    if(ex) then
      open(newunit=mf,file=trim(bqtxt))    
    else
      write(*,*)"[ERR] Missing mesh file"
      stop
    endif

    !!-------------------------Mesh Type 0--------------------!!
    write(9,'(" [MSG] meshRead Unit = ",I10)')mf
    write(9,'(" [MSG] C_K1, C_K2 = ",2I10)')C_K1,C_K2
    read(mf,*,end=11,err=11)bqtxt
    read(mf,*,end=11,err=11)b%nele,b%npl,b%nedg
    read(mf,*,end=11,err=11)bqtxt
    read(mf,*,end=11,err=11)b%nbnd,b%nbndtyp
    write(9,'(3a15)')"Elements","Linear Nodes","Edges"
    write(9,'(3i15)')b%nele,b%npl,b%nedg
    write(9,'(2a15)')"Bnd","BndTypes"
    write(9,'(2i15)')b%nbnd,b%nbndtyp

    !Assumption regarding number of quad Nodes
    b%npq=b%nedg
    b%npt=b%npl+b%npq

    !Nodes Read
    allocate(b%cor(b%npt,2))
    b%cor=-999  
    read(mf,*,end=11,err=11)bqtxt
    do i=1,b%npl
      read(mf,*,end=11,err=11)b%cor(i,1),b%cor(i,2)
    enddo
    write(9,*)"[MSG] Nodes read done"
    !Debug Comments
    ! do i=1,npoinl
    !   write(9,*)i,":",(coord(i,j),j=1,2),depth(i)
    ! enddo

    !Elements Read
    allocate(b%conn(b%nele,6))
    b%conn=0
    read(mf,*,end=11,err=11)bqtxt
    do i=1,b%nele
      read(mf,*,end=11,err=11)b%conn(i,3),&
        b%conn(i,1),b%conn(i,2)    
    enddo
    write(9,*)"[MSG] Elements read done"
    !Debug Comments
    ! do i=1,nelem
    !   write(9,*)i,":",(conn(i,j),j=1,6)
    ! enddo

    !Boundaries Read
    allocate(b%mabnd(b%nbnd,6))
    b%mabnd=0
    k=0;
    do j=1,b%nbndtyp
      read(mf,*,end=11,err=11)bqtxt
      read(mf,*,end=11,err=11)tmpi1,tmpi2  
      do i=k+1,k+tmpi2
        read(mf,*,end=11,err=11)b%mabnd(i,1:3)
        b%mabnd(i,4)=tmpi1      
      enddo
      k=k+tmpi2    
    enddo
    if(k.ne.b%nbnd) goto 14
    write(9,*) "[MSG] Boundaries read done" 
    ! ! Debug Comments 
    ! do i=1,nbnd
    !   write(9,*)i,":",(mabnd(i,j),j=1,6)
    ! enddo

    !Depth Read
    allocate(b%dep(b%npt))
    b%dep=-999
    read(mf,*,end=11,err=11)bqtxt
    do i=1,b%npl
      read(mf,*,end=11,err=11)b%dep(i)
    enddo
    write(9,*) "[MSG] Depth read done"
    !Debug Comments
    ! do i=1,npoinl
    !   write(9,*)i,depth(i)
    ! enddo


    goto 12
    !!-----------------------End Mesh Type 0------------------!!

    11 write(9,*) "[ERR] Check mesh file format"
    stop
    13 write(9,*) "[ERR] hex2dec error"
    stop
    14 write(9,*) "[ERR] Number of boundaries mismatch"
    stop
    12 close(mf)
    write(9,*)"[MSG] Done meshRead "
    write(9,*)
  end subroutine meshRead
!!---------------------------End meshRead--------------------------!!



!!-----------------------------femInit-----------------------------!!
  subroutine femInit(b)
  implicit none

    class(bsnqCase),intent(inout)::b
    integer(kind=C_k1)::nbndpoi
    integer(kind=C_k1),allocatable::tempia(:,:)
    
    b%maxNePoi=maxNePoi
    b%maxNeEle=maxNeEle
    
    !!-----------------Generate Middle Points-----------------!!
    allocate(b%p2p(b%npt,b%maxNePoi),b%p2e(b%npt,b%maxNeEle))
    allocate(b%npoisur(b%npt,3))
    b%p2p=0
    b%npoisur=0

    call middleNode(b%npl,b%npq,b%npt,b%nele,b%nedg,&
      b%maxNePoi,b%cor,b%dep,b%conn,b%p2p,b%npoisur)
    if(b%npq.ne.b%nedg) then
      write(9,*) "[ERR] Initial npoint assumption insufficiant"
      stop
    endif
    b%p2p=0
    b%npoisur=0

    !Debug Comments
    ! do i=1,npoint
    !   write(9,*)i,":",(coord(i,j),j=1,2),depth(i)
    ! enddo
    ! do i=1,nelem
    !   write(9,*)i,":",(conn(i,j),j=1,6)
    ! enddo
    write(9,*)"[MSG] Middle points generation done"
    write(9,'(" [INF] ",4A12)')"LinNode","QuadNode","TotNode",&
      "nEdges"
    write(9,'(" [---] ",4I12)')b%npl,b%npq,b%npt,b%nedg

    !Boundary sides middle nodes update
    !middle nodes stored in mabnd(i,5)
    call bndMiddleNodes(b%npl,b%npt,b%nele,b%nbnd,&
      b%conn,b%mabnd)
    write(9,*) "[MSG] Boundary sides middle point done"  
    !Debug Comments
    ! do i=1,nbnd
    !   write(9,*)i,":",(mabnd(i,j),j=1,6)
    ! enddo

    !Storing boundary nodes in bndNode
    nbndpoi=0
    allocate(tempia(b%npt,2))
    tempia=0
    do i=1,b%nbnd
      nq(1)=b%mabnd(i,1)
      nq(2)=b%mabnd(i,2)
      nq(3)=b%mabnd(i,5)
      tmpi4=b%mabnd(i,4)
      do j=1,3
        do k=1,nbndpoi
          if(tempia(k,1).eq.nq(j)) goto 31
        enddo
        nbndpoi=nbndpoi+1
        tempia(nbndpoi,1)=nq(j)
        tempia(nbndpoi,2)=tmpi4
        31 continue
        !Preferance for bndtype 14 over 12 or 13
        if((tempia(k,2).ne.14).and.(tmpi4.eq.14)) then
          tempia(k,2)=14
        end if
        !Preferance for bndtype 11 - higher
        if((tempia(k,2).ne.11).and.(tmpi4.eq.11)) then
          tempia(k,2)=11
        end if
        ! if((tempia(k,2).eq.12).and.(tmpi4.eq.13)) then
        !   tempia(k,2)=13
        ! end if
      enddo
    enddo
    allocate(b%bndP(nbndpoi),b%bndPT(nbndpoi))
    b%bndP=tempia(1:nbndpoi,1)
    b%bndPT=tempia(1:nbndpoi,2)
    b%nbndp=nbndpoi
    deallocate(tempia)
    allocate(tempia(nbndpoi,2))
    tempia(:,1)=b%bndP
    tempia(:,2)=b%bndPT
    b%bndPT=0
    call mergeSort(nbndpoi,nbndpoi,b%bndP)
    do i=1,nbndpoi
      tmpi1=b%bndP(i)
      do j=1,nbndpoi
        if(tmpi1.eq.tempia(j,1)) goto 35
      enddo
      35 b%bndPT(i)=tempia(j,2)
    enddo
    deallocate(tempia)
    ! write(9,*)nbndpoi
    ! do i=1,nbndpoi
    !   k=bndNode(i)
    !   write(9,*)i,":",bndNodeType(i),coord(k,1),coord(k,2)
    ! enddo
  !!-----------------End Generate Middle Points---------------!!

  !!----------------------Node Connectivity-------------------!!
  call nodeConnVSR(b%npt,b%nele,b%maxNePoi,b%maxNeEle,&
    b%conn,b%p2p,b%p2e,b%npoisur)
  do i=1,b%npt
    call mergeSort(b%maxNePoi,b%npoisur(i,1),b%p2p(i,:))
  enddo

  !Finding number of linear element nbhs
  do i=1,b%npt
    do j=1,b%npoisur(i,1)
      if(b%p2p(i,j).gt.b%npl) exit
      b%npoisur(i,3)=b%npoisur(i,3)+1
    enddo
  enddo

  !Debug Comments
  ! do i=1,npoint
  !   write(9,*) i,":",(poi2poi(i,j),j=1,maxNePoi)
  ! enddo
  ! do i=1,npoint
  !   write(9,*) i,":",(poi2ele(i,j),j=1,maxNeEle)
  ! enddo
  ! do i=1,npoint
  !   write(9,*) i,":",(npoisur(i,j),j=1,3)
  ! enddo
  write(9,*)"[MSG] Node Connectivity Done"
  !!--------------------End Node Connectivity-----------------!!

  !!--------------------------VSR Matrices--------------------!!
  !Storage allocations
  allocate(b%ivl(0:b%npt),b%ivq(0:b%npt))  
  b%ivl(0)=maxval(b%npoisur(:,3))+1
  b%ivq(0)=maxval(b%npoisur(:,1))+1
  allocate(b%linkl(b%ivl(0)*b%npt),b%linkq(b%ivq(0)*b%npt))
  b%linkl=0
  b%linkq=0

  !IV matrix linear and quadratic
  do i=1,b%npl
    b%ivl(i)=b%npoisur(i,3)+1
  enddo
  do i=b%npl+1,b%npt
    b%ivl(i)=b%npoisur(i,3)
  enddo
  ! write(9,*) ivl

  do i=1,b%npt
    b%ivq(i)=b%npoisur(i,1)+1
  enddo
  ! write(9,*) ivq

  !link matrix linear
  do i=1,b%npl
    k=(i-1)*b%ivl(0)
    do j=1,b%npoisur(i,3)
      b%linkl(k+j)=b%p2p(i,j)
    enddo
    b%linkl(k+b%ivl(i))=i
  enddo
  do i=b%npl+1,b%npt
    k=(i-1)*b%ivl(0)
    do j=1,b%npoisur(i,3)
      b%linkl(k+j)=b%p2p(i,j)
    enddo    
  enddo
  ! do i=1,npoint
  !   k=(i-1)*ivl(0)
  !   write(9,*) i,":",(linkl(k+j),j=1,ivl(i))
  ! enddo

  !link matrix quadratic
  do i=1,b%npt
    k=(i-1)*b%ivq(0)
    do j=1,b%npoisur(i,1)
      b%linkq(k+j)=b%p2p(i,j)
    enddo
    b%linkq(k+b%ivq(i))=i
  enddo
  ! write(9,*)linkq

  ! Full Matrices ivFull and linkf
  allocate(b%ivf(0:2*b%npt))
  b%ivf(0)=2*b%ivq(0)
  allocate(b%linkf(b%ivf(0)*2*b%npt))
  do i=1,b%npt
    b%ivf(i)=2*b%ivq(i)
    b%ivf(b%npt+i)=2*b%ivq(i)
  enddo
  
  do i=1,b%npt
    k=(i-1)*b%ivf(0)
    b%linkf(k+1:k+b%npoisur(i,1))=b%p2p(i,1:b%npoisur(i,1))
    b%linkf(k+b%npoisur(i,1)+1:k+(2*b%npoisur(i,1)))=b%npt &
      + b%p2p(i,1:b%npoisur(i,1))
    b%linkf(k+b%ivf(i)-1)=b%npt+i
    b%linkf(k+b%ivf(i))=i

    k=(i+b%npt-1)*b%ivf(0)
    b%linkf(k+1:k+b%npoisur(i,1))=b%p2p(i,1:b%npoisur(i,1))
    b%linkf(k+b%npoisur(i,1)+1:k+(2*b%npoisur(i,1)))=b%npt &
      + b%p2p(i,1:b%npoisur(i,1))
    b%linkf(k+b%ivf(i)-1)=i
    b%linkf(k+b%ivf(i))=b%npt+i
  enddo

  write(9,*)"[MSG] VSR storage matrices done"
  !!-----------------------End VSR Matrices-------------------!!

  !!---------------------Jacobian and Normals-----------------!!
  ! ShapeFnc Gauss-points JacobianQuad
  allocate(b%invJ(b%nele,5),b%bndS(b%nbnd,3))
  call jacbInvLin(b%npt,b%nele,b%conn,b%cor,b%invJ)      
  ! Priniting area using linear and quad jacb
  ! do iel=1,nelem
  !   tmpr1=invJ(iel,5)/2d0
  !   tmpr2=0d0
  !   do i=1,nGP
  !     tmpr2=tmpr2+gpW(i)*jacb(iel)%D(i)
  !   enddo
  !   !if(abs(tmpr2-tmpr1).lt.1e-10)cycle
  !   write(*,'(i10,3f15.8)')iel,tmpr1,tmpr2,abs(tmpr2-tmpr1)    
  ! enddo  
  ! stop

  ! Boundary side normals
  call bndSideInfo(b%npl,b%npt,b%nele,b%nbnd,b%cor,b%mabnd,&
    b%bndS)
  ! do i=1,nbnd
  !   write(9,*)i,":",bndSide(i,:),mabnd(i,1:5)
  ! enddo  

  !Boundary Node Normal
  allocate(b%bndPN(b%npt,2))
  call bndNodeNormal(b%npl,b%npt,b%nele,b%nbnd,b%cor,b%mabnd,&
    b%bndS,b%bndPN)
  ! do i=1,npoint
  !   write(9,*)i,":",bndNormal(i,1),bndNormal(i,2)
  ! enddo
  write(9,*)"[MSG] Done femInit"
  write(9,*)
  !!----------------End Jacobian and Normals------------------!!

  end subroutine femInit
!!---------------------------End femInit---------------------------!!



!!-----------------------------setRun------------------------------!!
  subroutine setRun(b)
  implicit none

    class(bsnqCase),intent(inout)::b
    integer(kind=C_K1)::mf

    !Input file open  
    bqtxt=trim(b%probname)//'.inp'
    inquire(file=trim(bqtxt),exist=ex)
    if(ex) then
      open(newunit=mf,file=trim(bqtxt))
    else
      write(9,*)"[ERR] Missing input file"
      stop
    endif

    write(9,'(" [MSG] setRun Unit = ",I10)')mf
    read(mf,*,end=81,err=81)bqtxt
    read(mf,*,end=81,err=81)b%resume
    read(mf,*,end=81,err=81)bqtxt
    read(mf,*,end=81,err=81)b%resumeFile  
    read(mf,*,end=81,err=81)bqtxt
    read(mf,*,end=81,err=81)b%dt     
    read(mf,*,end=81,err=81)bqtxt
    read(mf,*,end=81,err=81)b%nTSteps  
    read(mf,*,end=81,err=81)bqtxt
    read(mf,*,end=81,err=81)b%errLim
    read(mf,*,end=81,err=81)bqtxt
    read(mf,*,end=81,err=81)b%maxIter
    read(mf,*,end=81,err=81)bqtxt
    read(mf,*,end=81,err=81)i
    read(mf,*,end=81,err=81)bqtxt
    read(mf,*,end=81,err=81)tmpr1
    b%fileOut=int(tmpr1/b%dt,4)
    read(mf,*,end=81,err=81)bqtxt
    read(mf,*,end=81,err=81)tmpr1
    b%resumeOut=int(tmpr1/b%dt,4)  
    read(mf,*,end=81,err=81)bqtxt
    read(mf,*,end=81,err=81)b%presOn      
    read(mf,*,end=81,err=81)bqtxt
    read(mf,*,end=81,err=81)b%nthrd
    !read(mf,*,end=81,err=81)bqtxt
    !read(mf,*,end=81,err=81)tmpi1
    ! if(tmpi1.gt.0) then
    !   allocate(probe(0:tmpi1))
    !   probe(0)=tmpi1
    !   read(mf,*,end=81,err=81)probe(1:probe(0))
    ! else
    !   allocate(probe(0:1))
    !   probe(0)=0
    ! endif

    goto 82
    81 write(9,*) "[ERR] Check input file format"
    stop
    82 close(mf)

    write(9,*)"[MSG] Done setRun"
    write(9,*)
    
  end subroutine setRun
!!---------------------------End setRun----------------------------!!



!!-----------------------------initMat-----------------------------!!
  subroutine initMat(b)
  implicit none

    class(bsnqCase),intent(inout)::b

    i=b%npl
    j=b%npt
    i1=b%ivl(0)
    j1=b%ivq(0)

    allocate(b%pt0(j),b%qt0(j),b%et0(i),b%tDt0(j))
    allocate(b%pt1(j),b%qt1(j),b%et1(i),b%tDt1(j))
    allocate(b%pr(j),b%qr(j),b%er(i),b%tDr(j))
    allocate(b%wr(i),b%por(j))
    allocate(b%ur(j),b%vr(j),b%pbpr(j),b%qbpr(j))

    allocate(b%massW(i1*i),b%massE(i1*i))
    allocate(b%massP(j1*j),b%massQ(j1*j))    
    allocate(b%gBs1(j1*j),b%gBs2(j1*j))
    allocate(b%gBs3(j1*j),b%gBs4(j1*j))
    allocate(b%gCx(j1*i),b%gCy(j1*i),b%gDMat(i1*i))
    allocate(b%gBs5(i1*j),b%gBs6(i1*j))
    allocate(b%gGx(i1*j),b%gGy(i1*j),b%gNAdv(j1*j))
    allocate(b%gFBs1(j1*j),b%gFBs2(j1*j))
    allocate(b%gFBs3(j1*j),b%gFBs4(j1*j))

    b%por=1d0

    b%rTime=0d0
    b%et0=0d0
    b%pt0=0d0
    b%qt0=0d0

    call system_clock(b%sysClk(1),b%sysClk(2))
    b%sysRate=real(b%sysClk(2),C_K2)

    write(9,*)"[MSG] Done initMat"
    write(9,*)

  end subroutine initMat
!!---------------------------End initMat---------------------------!!



!!---------------------------statMatrices--------------------------!!
  subroutine statMatrices(b)
  implicit none

    class(bsnqCase),intent(inout)::b    

    call matrixSet1(b%npl,b%npt,b%nele,b%conn,b%ivl,b%ivq,&
      b%linkl,b%linkq,b%invJ,b%dep,b%por,b%massP,b%massE,&
      b%gBs1,b%gBs2,b%gBs3,b%gBs4,b%gCx,b%gCy,b%gDMat,&
      b%gBs5,b%gBs6)
    write(9,*)"[MSG] Done matrixSet1"

    call bndIntegral1(b%npl,b%npt,b%nele,b%nbnd,b%conn,b%mabnd,&
      b%ivl,b%ivq,b%linkl,b%linkq,b%invJ,b%bndS,b%dep,&
      b%gFBs1,b%gFBs2,b%gFBs3,b%gFBs4)
    write(9,*)"[MSG] Done bndIntegral1"

    b%massW=b%massE
    b%massQ=b%massP
    b%gBs1=b%gBs1+b%gFBs1
    b%gBs2=b%gBs2+b%gFBs2
    b%gBs3=b%gBs3+b%gFBs3
    b%gBs4=b%gBs4+b%gFBs4

    call b%destructR1

    write(9,*)"[MSG] Done statMatrices"
    write(9,*)

  end subroutine statMatrices
!!-------------------------End statMatrices------------------------!!



!!---------------------------dynaMatrices--------------------------!!
  subroutine dynaMatrices(b)
  implicit none

    class(bsnqCase),intent(inout)::b    

    call matrixSet2(b%npl,b%npt,b%nele,b%conn,b%ivl,b%ivq,&
      b%linkl,b%linkq,b%invJ,b%dep,b%por,b%tDr,b%ur,b%vr,&
      b%gGx,b%gGy,b%gNAdv)

    write(9,*)"[MSG] Done dynaMatrices"
    write(9,*)

  end subroutine dynaMatrices
!!-------------------------End dynaMatrices------------------------!!



!!----------------------------destructR1---------------------------!!
  subroutine destructR1(b)
  implicit none

    class(bsnqCase),intent(inout)::b    

    deallocate(b%gFBs1,b%gFBs2,b%gFBs3,b%gFBs4)

  end subroutine destructR1
!!--------------------------End destructR1-------------------------!!
end module bsnqModule