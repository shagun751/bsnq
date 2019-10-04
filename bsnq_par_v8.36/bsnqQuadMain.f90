!!---------------------- Version 8.x.x -----------------------!!
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

include 'subroutines/mods.f90'
include 'subroutines/bndIntegral.f90'
include 'subroutines/boundaryDifferential.f90'
include 'subroutines/geometry.f90'
include 'subroutines/gqMatrixSet1.f90'
include 'subroutines/matrixSet1.f90'
include 'subroutines/matrixSet2.f90'
include 'subroutines/mergeSort.f90'
include 'subroutines/nodeConnAll.f90'
include 'subroutines/outputN.f90'
include 'subroutines/porMatrices.f90'
include 'subroutines/resumeFile.f90'
include 'subroutines/shipFncs.f90'
include 'subroutines/waveCalculator_v2.f90'


program boussinesqQuad

use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE
use, intrinsic :: ISO_C_BINDING, only : C_CHAR, C_NULL_CHAR, C_LOC
use basicVars
use shapeFnc

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
    
  integer(kind=C_K1)::npoinl,npoinq,npoint,nelem,nbnd,nbndtyp,nedge    
  integer(kind=C_K1),allocatable::conn(:,:),mabnd(:,:),tempia(:,:)  
  integer(kind=C_K1),allocatable::poi2poi(:,:),poi2ele(:,:)
  integer(kind=C_K1),allocatable::npoisur(:,:),bndNode(:)
  integer(kind=C_K1),allocatable::bndNodeType(:),probe(:)
  integer(kind=C_K1),allocatable::outAbsN(:),outAbsP(:)
  integer(kind=C_K1)::nbndpoi,fileOut,resumeOut
  integer(kind=C_K1)::numTSteps,tStep,ier,maxIter,iter,gsIter
  integer(kind=C_K1)::driMeth,nthrd,nnzCSR1,nnzCSR2

  integer(kind=C_K1),allocatable::ivl(:),linkl(:),ivq(:),linkq(:)
  integer(kind=C_K1),allocatable::ivFull(:),linkf(:)
  integer(kind=C_INT), allocatable, target :: ivCSR1(:),jvCSR1(:)
  integer(kind=C_INT), allocatable, target :: ivCSR2(:),jvCSR2(:)

  real(kind=C_K2),allocatable::coord(:,:),depth(:),totDe(:)
  real(kind=C_K2),allocatable::invJ(:,:),mass1(:),mass2(:)
  real(kind=C_K2),allocatable::gBs1(:),gBs2(:),gBs3(:),gBs4(:)
  real(kind=C_K2),allocatable::gM1Bs1(:),gM1Bs4(:),mass2eta(:)
  real(kind=C_K2),allocatable::gCxF(:),gCyF(:),mass2A(:)
  real(kind=C_K2),allocatable::gAx(:),gAy(:),gDMat(:),gdNdn(:)
  real(kind=C_K2),allocatable::p(:),q(:),u(:),v(:),eta(:)
  real(kind=C_K2),allocatable::pt1(:),qt1(:),etat1(:),wt1(:)
  real(kind=C_K2),allocatable::pdt(:),qdt(:),wt2(:)
  real(kind=C_K2),allocatable::pdtt1(:),qdtt1(:),etadtt1(:)
  real(kind=C_K2),allocatable::pdtt2(:),qdtt2(:),etadtt2(:)
  real(kind=C_K2),allocatable::gNAdv(:),gGxElev(:),gGyElev(:)
  real(kind=C_K2),allocatable::bndSide(:,:),bndNormal(:,:)
  real(kind=C_K2),allocatable::gFp(:),gFq(:),gFn(:),gFw(:)
  real(kind=C_K2),allocatable::g1(:),g2(:),gRun(:)
  real(kind=C_K2),allocatable::aFull(:)
  real(kind=C_K2),allocatable::pRun(:),qRun(:),etaRun(:)
  real(kind=C_K2),allocatable::etadtRun(:)
  real(kind=C_K2),allocatable::etaMin(:),etaMax(:),por(:)
  real(kind=C_K2),allocatable::porLam(:),porTurX(:),porTurY(:)
  real(kind=C_K2),allocatable::toux(:),touy(:),touMass(:)
  real(kind=C_K2),allocatable::velAbs(:,:)
  real(kind=C_K2),allocatable::etaAbs(:)
  real(kind=C_K2),allocatable::presDiv(:,:),presK(:)
  real(kind=C_K2),allocatable::rowMaxA(:),rowMaxB(:),rowMaxC(:)
  real(kind=C_K2),allocatable::outAbsV(:,:),etaAbsC(:),velAbsC(:)
  real(kind=C_K2),allocatable::tempra(:,:),porH(:),porV(:)
  real(kind=C_K2)::tempr2(5),cnst(10),errLim,errOut
  real(kind=C_K2)::alpha,alpha2,rTime,dt
  real(kind=C_K2)::waveT(5),waveL(5),waveK(5),waveW(5),waveA(5)
  real(kind=C_K2)::stoneCnst(10),waveInfo(5,5)
  real(kind=C_K2)::outAbs(10),resnorm  

  type(jacbType),allocatable::jacb(:)
  
  real(kind=C_DOUBLE), allocatable, target :: aFullCSR(:),mass2CSR(:),mass2etaCSR(:)
  real(kind=C_DOUBLE), allocatable, target :: g3(:),g4(:),bFull(:)
  real(kind=C_DOUBLE), allocatable, target :: w(:),xFull(:),etadt(:)  

  character(len=100)::probname,text,resumeFile
  
  logical::ex,resume,presOn,outAbsOn,porOn  
  logical,allocatable::templa(:)  

!!--------------------------All Constants---------------------!!
  !!---------------------Initialisations-----------------!!  
  nbnd=0  

  !!------------------------Constants--------------------!!
  cnst(1)=grav                    !gravity
  cnst(2)=BsqC                    !Bsnq constant
  cnst(3)=BsqC+(1d0/3d0)          !dervied 1
  cnst(4)=2d0*BsqC+(1d0/3d0)      !derived 2
  cnst(5)=2d0*BsqC+(1d0/2d0)      !derived 3
  cnst(6)=grav*BsqC               !derived 4
  cnst(7)=pi                      !pi
  cnst(8)=rhoW                    !Density of water
  cnst(9)=0.015d0                 !friction factor = 1/2*rho*Cb

  !!---------------------Stone Constants-----------------!!
  stoneCnst(1)=1100d0             !alpha0 780 to 1500
  stoneCnst(2)=1.50d0             !beta0  1.8 to 3.6
  stoneCnst(3)=1e-6               !kinematic viscosity of water
  stoneCnst(4)=0.02d0             !Stone size
  stoneCnst(5)=stoneCnst(1)*stoneCnst(3)/(stoneCnst(4)**2)
  stoneCnst(6)=stoneCnst(2)/stoneCnst(4)


  call getarg(1,probname)  
  if(len_trim(probname).lt.1) then
    write(*,'(A)',advance='no')"Enter Problem Name:"
    read(*,*)probname
  endif
  write(*,*)"Problem Name: "//probname(1:len_trim(probname))
!!------------------------End All Constants-------------------!!

!!-------------------------Mesh File Reading------------------!!
  !Opening mesh file  
  text=probname(1:len_trim(probname))//'.plt'
  inquire(file=text(1:len_trim(text)),exist=ex)
  if(ex) then
    open(mafi(1),file=text(1:len_trim(text)))    
  else
    write(*,*)"[Error] Missing mesh file"
    stop
  endif

  text=probname(1:len_trim(probname))//'.rout'
  open(9,file=text(1:len_trim(text)))

  !!-------------------------Mesh Type 0--------------------!! 
  read(mafi(1),*,end=11,err=11)text
  read(mafi(1),*,end=11,err=11)nelem,npoinl,nedge
  read(mafi(1),*,end=11,err=11)text
  read(mafi(1),*,end=11,err=11)nbnd,nbndtyp
  write(9,'(3a15)')"Elements","Linear Nodes","Edges"
  write(9,'(3i15)')nelem,npoinl,nedge
  write(9,'(2a15)')"Bnd","BndTypes"
  write(9,'(2i15)')nbnd,nbndtyp

  !Assumption regarding number of quad Nodes
  npoinq=nedge
  npoint=npoinl+npoinq

  !Nodes Read
  allocate(coord(npoint,2))
  coord=-999  
  read(mafi(1),*,end=11,err=11)text
  do i=1,npoinl
    read(mafi(1),*,end=11,err=11)coord(i,1),coord(i,2)
  enddo
  write(9,*)"Nodes read done"
  !Debug Comments
  ! do i=1,npoinl
  !   write(9,*)i,":",(coord(i,j),j=1,2),depth(i)
  ! enddo

  !Elements Read
  allocate(conn(nelem,6))
  conn=0
  read(mafi(1),*,end=11,err=11)text
  do i=1,nelem    
    read(mafi(1),*,end=11,err=11)conn(i,3),conn(i,1),conn(i,2)    
  enddo
  write(9,*)"Elements read done"
  !Debug Comments
  ! do i=1,nelem
  !   write(9,*)i,":",(conn(i,j),j=1,6)
  ! enddo

  !Boundaries Read
  allocate(mabnd(nbnd,6))
  mabnd=0
  k=0;
  do j=1,nbndtyp
    read(mafi(1),*,end=11,err=11)text
    read(mafi(1),*,end=11,err=11)tmpi1,tmpi2  
    do i=k+1,k+tmpi2
      read(mafi(1),*,end=11,err=11)mabnd(i,1:3)
      mabnd(i,4)=tmpi1      
    enddo
    k=k+tmpi2    
  enddo
  if(k.ne.nbnd) goto 14
  write(9,*) "Boundaries read done" 
  ! ! Debug Comments 
  ! do i=1,nbnd
  !   write(9,*)i,":",(mabnd(i,j),j=1,6)
  ! enddo

  !Depth Read
  allocate(depth(npoint))
  depth=-999
  read(mafi(1),*,end=11,err=11)text
  do i=1,npoinl
    read(mafi(1),*,end=11,err=11)depth(i)
  enddo
  write(9,*) "Depth read done"
  !Debug Comments
  ! do i=1,npoinl
  !   write(9,*)i,depth(i)
  ! enddo

  goto 12
  !!-----------------------End Mesh Type 0------------------!!

  11 write(9,*) "[Error] Check mesh file format"
  stop
  13 write(9,*) "[Error] hex2dec error"
  stop
  14 write(9,*) "[Error] Number of boundaries mismatch"
  stop
  12 close(mafi(1))
  write(9,*) "Mesh file read successful"
!!---------------------End Mesh File Reading------------------!!

!!--------------------Generate Middle Points------------------!!
  allocate(poi2poi(npoint,maxNePoi),poi2ele(npoint,maxNeEle))
  allocate(npoisur(npoint,3))
  poi2poi=0
  npoisur=0

  call middleNode(npoinl,npoinq,npoint,nelem,nedge,maxNePoi,&
    coord,depth,conn,poi2poi,npoisur)
  if(npoinq.ne.nedge) then
    write(9,*) "[Error] Initial npoint assumption insufficiant"
    stop
  endif
  poi2poi=0
  npoisur=0

  !Debug Comments
  ! do i=1,npoint
  !   write(9,*)i,":",(coord(i,j),j=1,2),depth(i)
  ! enddo
  ! do i=1,nelem
  !   write(9,*)i,":",(conn(i,j),j=1,6)
  ! enddo
  write(9,*)"Middle points generation done"
  write(9,*)"LinNode, QuadNode, TotNode, nEdges"
  write(9,*) npoinl,npoinq,npoint,nedge

  !Boundary sides middle nodes update
  !middle nodes stored in mabnd(i,5)
  call bndMiddleNodes(npoinl,npoint,nelem,nbnd,conn,mabnd)
  write(9,*) "Boundary sides middle point done"  
  !Debug Comments
  ! do i=1,nbnd
  !   write(9,*)i,":",(mabnd(i,j),j=1,6)
  ! enddo

  !Storing boundary nodes in bndNode
  nbndpoi=0
  allocate(tempia(npoint,2))
  tempia=0
  do i=1,nbnd
    nq(1)=mabnd(i,1)
    nq(2)=mabnd(i,2)
    nq(3)=mabnd(i,5)
    tmpi4=mabnd(i,4)
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
  allocate(bndNode(nbndpoi),bndNodeType(nbndpoi))
  bndNode=tempia(1:nbndpoi,1)
  bndNodeType=tempia(1:nbndpoi,2)
  deallocate(tempia)
  allocate(tempia(nbndpoi,2))
  tempia(:,1)=bndNode
  tempia(:,2)=bndNodeType
  bndNodeType=0
  call mergeSort(nbndpoi,nbndpoi,bndNode)
  do i=1,nbndpoi
    tmpi1=bndNode(i)
    do j=1,nbndpoi
      if(tmpi1.eq.tempia(j,1)) goto 35
    enddo
    35 bndNodeType(i)=tempia(j,2)
  enddo
  deallocate(tempia)
  ! write(9,*)nbndpoi
  ! do i=1,nbndpoi
  !   k=bndNode(i)
  !   write(9,*)i,":",bndNodeType(i),coord(k,1),coord(k,2)
  ! enddo
!!------------------Done Generate Middle Points---------------!!

!!-----------------------Node Connectivity--------------------!!
  call nodeConnVSR(npoint,nelem,maxNePoi,maxNeEle,conn,&
    poi2poi,poi2ele,npoisur)
  do i=1,npoint
    call mergeSort(maxNePoi,npoisur(i,1),poi2poi(i,:))
  enddo

  !Finding number of linear element nbhs
  do i=1,npoint
    do j=1,npoisur(i,1)
      if(poi2poi(i,j).gt.npoinl) exit
      npoisur(i,3)=npoisur(i,3)+1
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
  write(9,*)"Node Connectivity Done"
!!---------------------Done Node Connectivity-----------------!!

!!---------------------------VSR Matrices---------------------!!
  !Storage allocations
  allocate(ivl(0:npoint),ivq(0:npoint))  
  ivl(0)=maxval(npoisur(:,3))+1
  ivq(0)=maxval(npoisur(:,1))+1
  allocate(linkl(ivl(0)*npoint),linkq(ivq(0)*npoint))
  linkl=0
  linkq=0

  !IV matrix linear and quadratic
  do i=1,npoinl
    ivl(i)=npoisur(i,3)+1
  enddo
  do i=npoinl+1,npoint
    ivl(i)=npoisur(i,3)
  enddo
  ! write(9,*) ivl

  do i=1,npoint
    ivq(i)=npoisur(i,1)+1
  enddo
  ! write(9,*) ivq

  !link matrix linear
  do i=1,npoinl
    k=(i-1)*ivl(0)
    do j=1,npoisur(i,3)
      linkl(k+j)=poi2poi(i,j)
    enddo
    linkl(k+ivl(i))=i
  enddo
  do i=npoinl+1,npoint
    k=(i-1)*ivl(0)
    do j=1,npoisur(i,3)
      linkl(k+j)=poi2poi(i,j)
    enddo    
  enddo
  ! do i=1,npoint
  !   k=(i-1)*ivl(0)
  !   write(9,*) i,":",(linkl(k+j),j=1,ivl(i))
  ! enddo

  !link matrix quadratic
  do  i=1,npoint
    k=(i-1)*ivq(0)
    do j=1,npoisur(i,1)
      linkq(k+j)=poi2poi(i,j)
    enddo
    linkq(k+ivq(i))=i
  enddo
  ! write(9,*)linkq

  ! Full Matrices ivFull and linkf
  allocate(ivFull(0:2*npoint))
  ivFull(0)=2*ivq(0)
  allocate(linkf(ivFull(0)*2*npoint))
  do i=1,npoint
    ivFull(i)=2*ivq(i)
    ivFull(npoint+i)=2*ivq(i)
  enddo
  
  do i=1,npoint
    k=(i-1)*ivFull(0)
    linkf(k+1:k+npoisur(i,1))=poi2poi(i,1:npoisur(i,1))
    linkf(k+npoisur(i,1)+1:k+(2*npoisur(i,1)))=npoint &
      + poi2poi(i,1:npoisur(i,1))
    linkf(k+ivFull(i)-1)=npoint+i
    linkf(k+ivFull(i))=i

    k=(i+npoint-1)*ivFull(0)
    linkf(k+1:k+npoisur(i,1))=poi2poi(i,1:npoisur(i,1))
    linkf(k+npoisur(i,1)+1:k+(2*npoisur(i,1)))=npoint &
      + poi2poi(i,1:npoisur(i,1))
    linkf(k+ivFull(i)-1)=i
    linkf(k+ivFull(i))=npoint+i
  enddo

  write(9,*)"VSR storage matrices done"
!!------------------------Done VSR Matrices-------------------!!

!!-------------------------Jacobian and Normals---------------!!
  ! ShapeFnc Gauss-points JacobianQuad
  allocate(invJ(nelem,5),bndSide(nbnd,3),jacb(nelem))
  call jacbInvLin(npoint,nelem,conn,coord,invJ)    
  call initGaussPoi
  call initShapeFnc
  call jacbInvQuad(npoint,nelem,conn,coord,jacb)
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
  call bndSideInfo(npoinl,npoint,nelem,nbnd,coord,mabnd,&
    bndSide)
  ! do i=1,nbnd
  !   write(9,*)i,":",bndSide(i,:),mabnd(i,1:5)
  ! enddo
  write(9,*)"Boundary Geometry parameters done"

  !Boundary Node Normal
  allocate(bndNormal(npoint,2))
  call bndNodeNormal(npoinl,npoint,nelem,nbnd,coord,mabnd,&
    bndSide,bndNormal)
  ! do i=1,npoint
  !   write(9,*)i,":",bndNormal(i,1),bndNormal(i,2)
  ! enddo
!!-------------------Done Jacobian and Normals----------------!!

!!------------------------Allocations-------------------------!!
  allocate(mass1(ivq(0)*npoint),mass2(ivl(0)*npoinl))  
  allocate(gBs1(ivq(0)*npoint),gBs2(ivq(0)*npoint))
  allocate(gBs3(ivq(0)*npoint),gBs4(ivq(0)*npoint))
  allocate(gM1Bs1(ivq(0)*npoint),gM1Bs4(ivq(0)*npoint))
  allocate(gCxF(ivq(0)*npoinl),gCyF(ivq(0)*npoinl))
  allocate(gAx(ivl(0)*npoint),gAy(ivl(0)*npoint))
  allocate(gDMat(ivl(0)*npoinl),mass2eta(ivl(0)*npoinl))

  allocate(p(npoint),q(npoint),u(npoinl),v(npoinl))
  allocate(eta(npoinl),w(npoinl),totDe(npoinl))
  allocate(gFp(npoint),gFq(npoint),gFn(npoinl),gFw(npoinl))

  allocate(gNAdv(ivq(0)*npoint),mass2A(ivl(0)*npoinl))
  allocate(gdNdn(ivl(0)*npoinl))
  allocate(gGxElev(ivl(0)*npoint),gGyElev(ivl(0)*npoint))

  allocate(g1(npoint),g2(npoint),g3(npoinl),g4(npoinl))
  allocate(etadt(npoinl),pdt(npoint),qdt(npoint),xFull(2*npoint))
  allocate(pdtt1(npoint),qdtt1(npoint),gRun(npoint))
  allocate(etadtt1(npoinl),etadtt2(npoinl),wt1(npoinl))
  allocate(pdtt2(npoint),qdtt2(npoint),wt2(npoinl))
  allocate(etat1(npoinl),pt1(npoint),qt1(npoint))
  allocate(pRun(npoint),qRun(npoint),etaRun(npoinl))
  allocate(etadtRun(npoinl))
  allocate(etaMin(npoinl),etaMax(npoinl))
  allocate(aFull(ivFull(0)*2*npoint),bFull(2*npoint))
  allocate(por(npoinl),porLam(ivq(0)*npoint))
  allocate(porH(npoinl),porV(npoinl))
  allocate(porTurX(ivq(0)*npoint),porTurY(ivq(0)*npoint))
  allocate(toux(npoinl),touy(npoinl),touMass(npoint*ivl(0)))
  allocate(presDiv(npoint,2),velAbs(npoint,2),etaAbs(npoinl))
  allocate(presK(ivq(0)*npoint))
  allocate(etaAbsC(npoinl),velAbsC(npoint))
  allocate(rowMaxA(npoinl),rowMaxB(npoinl),rowMaxC(2*npoint))  
!!-----------------------End Allocations----------------------!!
    
!!-----------------------------Loop---------------------------!!
  !!----------------------Settings-----------------------!!
  !Input file open  
  text=probname(1:len_trim(probname))//'.inp'
  inquire(file=text(1:len_trim(text)),exist=ex)
  if(ex) then
    open(mafi(5),file=text(1:len_trim(text)))    
  else
    write(9,*)"[Error] Missing input file"
    stop
  endif

  read(mafi(5),*,end=81,err=81)text
  read(mafi(5),*,end=81,err=81)resume
  read(mafi(5),*,end=81,err=81)text
  read(mafi(5),*,end=81,err=81)resumeFile  
  read(mafi(5),*,end=81,err=81)text
  read(mafi(5),*,end=81,err=81)dt     
  read(mafi(5),*,end=81,err=81)text
  read(mafi(5),*,end=81,err=81)numTSteps  
  read(mafi(5),*,end=81,err=81)text
  read(mafi(5),*,end=81,err=81)errLim
  read(mafi(5),*,end=81,err=81)text
  read(mafi(5),*,end=81,err=81)maxIter
  read(mafi(5),*,end=81,err=81)text
  read(mafi(5),*,end=81,err=81)driMeth
  read(mafi(5),*,end=81,err=81)text
  read(mafi(5),*,end=81,err=81)tmpr1
  fileOut=int(tmpr1/dt,4)
  read(mafi(5),*,end=81,err=81)text
  read(mafi(5),*,end=81,err=81)tmpr1
  resumeOut=int(tmpr1/dt,4)  
  read(mafi(5),*,end=81,err=81)text
  read(mafi(5),*,end=81,err=81)presOn      
  read(mafi(5),*,end=81,err=81)text
  read(mafi(5),*,end=81,err=81)nthrd
  read(mafi(5),*,end=81,err=81)text
  read(mafi(5),*,end=81,err=81)tmpi1
  if(tmpi1.gt.0) then
    allocate(probe(0:tmpi1))
    probe(0)=tmpi1
    read(mafi(5),*,end=81,err=81)probe(1:probe(0))
  else
    allocate(probe(0:1))
    probe(0)=0
  endif
  
  !!--------------------End Settings---------------------!!

  !!----------------Variable Initialisation--------------!!
  rTime=0d0
  eta=0d0
  w=0d0
  p=0d0
  q=0d0
  etadt=0d0
  pdt=0d0
  qdt=0d0
  wt1=0d0
  etadtt1=0d0
  pdtt1=0d0
  qdtt1=0d0
  gFp=0d0
  gFq=0d0
  gFn=0d0
  gFw=0d0
  presDiv=0d0
  etaMax=eta
  etaMin=eta  
  call system_clock(sysClk(1),sysRaInt)
  sysRate=real(sysRaInt,8)
  !!--------------End Variable Initialisation------------!!

  !!--------------Inlet Boundary Condition---------------!!
  waveT=1;
  waveL=1;
  waveA=1;
  waveInfo=0d0;
  !!waveInfo
  !![CentreX CentreY WaveAngleRad CosAng SinAng]

  ! Inlet Wave Characteristics
  read(mafi(5),*,end=81,err=81)text
  read(mafi(5),*,end=81,err=81)waveT(1),waveA(1)
  read(mafi(5),*,end=81,err=81)text
  read(mafi(5),*,end=81,err=81)waveInfo(1,1:3)
  waveInfo(1,3)=waveInfo(1,3)*pi/180d0;

  !!Assuming all inlet nodes are at the same water depth
  tmpr1=-999d0
  do i=1,nbndpoi
    if(bndNodeType(i).eq.11) then
      tmpr1=depth(bndNode(i))
      exit
    endif
  enddo
  if(tmpr1.ne.-999d0) then
    do i=1,1
      call waveCalculator(grav,pi,errLim,waveT(i),&
        tmpr1,waveL(i))
      write(9,'(a17,4f14.4)')"[INF] L,T,d,C :",waveL(i),&
        waveT(i),tmpr1,waveL(i)/waveT(i)
      if(waveL(i).lt.0d0)then
        write(9,'(" [ERR] Check waveLength")')
      endif
    enddo
  endif
  waveW=2*pi/waveT;
  waveK=2*pi/waveL;
  waveInfo(:,4)=dcos(waveInfo(:,3))
  waveInfo(:,5)=dsin(waveInfo(:,3))

  ! OutletAbs Inititalisation
  !outAbsOn=.true.
  ! outAbs(1)=30d0  !X
  ! outAbs(2)=0d0   !Y
  ! outAbs(3)=6.576d0   !delta
  ! outAbs(4)=30d0/waveT(1)
  ! outAbs(5)=exp(1d0)-1d0
  ! outAbs(6)=135d0*pi/180d0   !angle
  ! outAbs(6)=tan(outAbs(6))  !m
  ! outAbs(7)=outAbs(2)-(outAbs(6)*outAbs(1)) !c
  !!------------End Inlet Boundary Condition-------------!!

  !!------------Solitary Wave Inititalisation------------!!
  ! ! ! Solitary wave fnc2
  ! tmpr1=0.45d0
  ! tmpr2=0.045d0
  ! !depth=tmpr1
  ! tmpr3=dsqrt(grav*(tmpr1+tmpr2))
  ! tmpr4=dsqrt(3*tmpr2/(4*(tmpr1**3)))
  ! do i=1,npoint
  !   if((coord(i,1).ge.3d0).and.(coord(i,1).le.19d0)) then
  !     tmpr5=tmpr2/(dcosh(tmpr4*(coord(i,1)-(11d0)))**2)
  !     p(i)=tmpr3*tmpr5
  !     if(i.le.npoinl) then
  !       eta(i)=tmpr5
  !     endif
  !   endif
  ! enddo
  !!----------End Solitary Wave Inititalisation----------!!

  !!----------------OutletAbs Coefficient----------------!!
  read(mafi(5),*,end=81,err=81)text
  read(mafi(5),*,end=81,err=81)text
  read(mafi(5),*,end=81,err=81)outAbsOn
  read(mafi(5),*,end=81,err=81)text
  read(mafi(5),*,end=81,err=81)tmpi1
  if(tmpi1.le.0) then 
    outAbsOn=.false.
  else 
    allocate(outAbsN(0:tmpi1),outAbsV(tmpi1,6))
    outAbsN(0)=tmpi1
    do i=1,outAbsN(0)
      read(mafi(5),*,end=81,err=81)outAbsN(i),outAbsV(i,1:2)    
    enddo
    outAbsV(:,3)=30d0/waveT(1)
    outAbsV(:,4)=exp(1d0)-1d0
  endif
  
  if(outAbsOn) then
    !Find Sponge layer points
    ! outAbsP matirx defninition
    ! outAbsP(-1:NAbsPoi)
    ! outAbsP(-1)  = NAbsPoi Lin only
    ! outAbsP(0)   = NAbdPoi All = NAbsPoi
    ! outAbsP(1:outAbsP(-1)) = sponge layer lin node indices
    ! outAbsP(outAbsP(-1)+1:outAbsP(0)) = sponge layer quad n
    tmpi1=0    
    tmpi2=0    
    allocate(tempia(npoint,1),templa(npoint))
    etaAbsC=0d0
    velAbsC=0d0
    templa=.true.

    ! OutletAbs Coeff (lin poi only)  East West
    do i2=1,outAbsN(0)
      select case (outAbsN(i2))
        case (2) !West
          do i=1,npoinl
            if(coord(i,1).lt.outAbsV(i2,1)) then
              tmpr1=coord(i,1)-outAbsV(i2,1)
              tmpr2=outAbsV(i2,3) &
                *(exp((tmpr1/outAbsV(i2,2))**2)-1d0) &
                /outAbsV(i2,4)            
              etaAbsC(i)=tmpr2              
              velAbsC(i)=tmpr2              

              templa(i)=.false.
              tmpi1=tmpi1+1
              tempia(tmpi1,1)=i
            endif
          enddo        

        case (4) !East
          do i=1,npoinl
            if(coord(i,1).gt.outAbsV(i2,1)) then
              tmpr1=coord(i,1)-outAbsV(i2,1)
              tmpr2=outAbsV(i2,3) &
                *(exp((tmpr1/outAbsV(i2,2))**2)-1d0) &
                /outAbsV(i2,4)            
              etaAbsC(i)=tmpr2              
              velAbsC(i)=tmpr2              

              templa(i)=.false.
              tmpi1=tmpi1+1
              tempia(tmpi1,1)=i
            endif
          enddo                      
      end select      
    enddo

    ! OutletAbs Coeff (lin poi only)  North South
    do i2=1,outAbsN(0)
      select case (outAbsN(i2))
        case (1) !North
          do i=1,npoinl
            if(coord(i,2).gt.outAbsV(i2,1)) then
              tmpr1=coord(i,2)-outAbsV(i2,1)
              tmpr2=outAbsV(i2,3) &
                *(exp((tmpr1/outAbsV(i2,2))**2)-1d0) &
                /outAbsV(i2,4)            

              if(templa(i))then
                etaAbsC(i)=tmpr2              
                velAbsC(i)=tmpr2              

                templa(i)=.false.
                tmpi1=tmpi1+1
                tempia(tmpi1,1)=i
              endif
            endif
          enddo        

        case (3) !South
          do i=1,npoinl
            if(coord(i,2).lt.outAbsV(i2,1)) then
              tmpr1=coord(i,2)-outAbsV(i2,1)
              tmpr2=outAbsV(i2,3) &
                *(exp((tmpr1/outAbsV(i2,2))**2)-1d0) &
                /outAbsV(i2,4)            
                
              if(templa(i))then
                etaAbsC(i)=tmpr2              
                velAbsC(i)=tmpr2              

                templa(i)=.false.
                tmpi1=tmpi1+1
                tempia(tmpi1,1)=i
              endif
            endif
          enddo        
      end select      
    enddo

    ! OutletAbs Coeff (quad poi only) East West
    tmpi2=tmpi1
    do i2=1,outAbsN(0)
      select case (outAbsN(i2))
        case (2)
          do i=npoinl+1,npoint
            if(coord(i,1).lt.outAbsV(i2,1)) then
              tmpr1=coord(i,1)-outAbsV(i2,1)
              tmpr2=outAbsV(i2,3) &
                *(exp((tmpr1/outAbsV(i2,2))**2)-1d0) &
                /outAbsV(i2,4)                          
              velAbsC(i)=tmpr2              

              templa(i)=.false.
              tmpi2=tmpi2+1
              tempia(tmpi2,1)=i
            endif
          enddo        

        case (4)
          do i=npoinl+1,npoint
            if(coord(i,1).gt.outAbsV(i2,1)) then
              tmpr1=coord(i,1)-outAbsV(i2,1)
              tmpr2=outAbsV(i2,3) &
                *(exp((tmpr1/outAbsV(i2,2))**2)-1d0) &
                /outAbsV(i2,4)                          
              velAbsC(i)=tmpr2              

              templa(i)=.false.
              tmpi2=tmpi2+1
              tempia(tmpi2,1)=i
            endif
          enddo                
      end select      
    enddo

    ! OutletAbs Coeff (quad poi only)  North South
    do i2=1,outAbsN(0)
      select case (outAbsN(i2))
        case (1) !North
          do i=npoinl+1,npoint
            if(coord(i,2).gt.outAbsV(i2,1)) then
              tmpr1=coord(i,2)-outAbsV(i2,1)
              tmpr2=outAbsV(i2,3) &
                *(exp((tmpr1/outAbsV(i2,2))**2)-1d0) &
                /outAbsV(i2,4)            

              if(templa(i))then                
                velAbsC(i)=tmpr2              

                templa(i)=.false.
                tmpi2=tmpi2+1
                tempia(tmpi2,1)=i
              endif
            endif
          enddo        

        case (3) !South
          do i=npoinl+1,npoint
            if(coord(i,2).lt.outAbsV(i2,1)) then
              tmpr1=coord(i,2)-outAbsV(i2,1)
              tmpr2=outAbsV(i2,3) &
                *(exp((tmpr1/outAbsV(i2,2))**2)-1d0) &
                /outAbsV(i2,4)            

              if(templa(i))then                
                velAbsC(i)=tmpr2              

                templa(i)=.false.
                tmpi2=tmpi2+1
                tempia(tmpi2,1)=i
              endif
            endif
          enddo        
      end select      
    enddo

    allocate(outAbsP(-1:tmpi2))
    outAbsP(-1)=tmpi1
    outAbsP(0)=tmpi2
    outAbsP(1:outAbsP(0))=tempia(1:tmpi2,1)
    deallocate(tempia,templa)
    ! write(9,*)outAbsP(-1:0)
    ! write(9,*)outAbsP(1:outAbsP(-1))
    ! write(9,*)outAbsP(outAbsP(-1)+1:outAbsP(0))
    write(9,*)'[Msg] Outlet Sponge layer initialised'
  endif  
  !!--------------End OutletAbs Coefficient--------------!!

  !!---------------Porosity Initialisation---------------!!
  ! Porosity
  read(mafi(5),*,end=81,err=81)text
  read(mafi(5),*,end=81,err=81)text
  read(mafi(5),*,end=81,err=81)porOn
  read(mafi(5),*,end=81,err=81)text
  read(mafi(5),*,end=81,err=81)tmpi1  

  porV=1.0d0
  porH=0d0
  por=1.0d0  

  if(tmpi1.lt.0) then
    porOn=.false.

  else if(tmpi1.eq.0) then
    if(porOn)then      
      text=probname(1:len_trim(probname))//'_por.plt'
      inquire(file=text(1:len_trim(text)),exist=ex)
      if(ex) then
        open(mafi(6),file=text(1:len_trim(text)))        
      else
        write(9,*)"[Err] Missing porosity value file"
        stop
      endif

      read(mafi(6),*,end=91,err=91)text
      do i=1,npoinl
        read(mafi(6),*,end=91,err=91)porV(i),porH(i)
      enddo

      goto 92
      91 write(9,*) "[Err] Check porosity file format"
      stop
      92 close(mafi(6))
      write(9,*)'[Msg] Porosity file read'  
    endif

  else    
    allocate(tempra(tmpi1,5))
    do i2=1,tmpi1
      read(mafi(5),*,end=81,err=81)tempra(i2,1:5)      
    enddo
    if(porOn) then
      do i=1,npoinl
        do i2=1,tmpi1
          if((coord(i,1).ge.tempra(i2,2)).and.(coord(i,1).le.tempra(i2,4))) then
            if((coord(i,2).ge.tempra(i2,3)).and.(coord(i,2).le.tempra(i2,5))) then
              porV(i)=tempra(i2,1)
            endif
          endif
        enddo
      enddo    
      write(9,*)'[Msg] Porosity Initialisation done'  
    endif
    deallocate(tempra)
  endif  
  
  ! tmpr4=0.5d0
  ! tmpr1=4.995d0
  ! tmpr2=5.305d0
  ! do i=1,npoinl
  !   if((coord(i,1).ge.tmpr1).and.(coord(i,1).le.tmpr2)) then
  !       por(i)=tmpr4
  !   endif
  ! enddo
  !!-------------End Porosity Initialisation-------------!!

  goto 82
  81 write(9,*) "[Error] Check input file format"
  stop
  82 close(mafi(5))

  ! Constant Matrices
  call massMatrices(npoint,npoinl,nelem,conn,ivl,ivq,linkl,linkq,&
    invJ,mass1,mass2,porLam,por,stoneCnst,touMass)
    
  call matrixSet1(npoinl,npoint,nelem,conn,ivl,ivq,linkl,linkq,&
    invJ,depth,gBs1,gBs2,gBs3,gBs4,gCxF,gCyF,gAx,gAy,&
    gDMat,por,cnst)

  call etaBC(npoinl,npoint,nelem,nbndpoi,maxNeEle,conn,poi2ele,&
    bndNode,ivl,linkl,invJ,bndNormal,gdNdn)

  gM1Bs1=gBs1+mass1
  gM1Bs4=gBs4+mass1
  mass2A=mass2

  if(driMeth.eq.2) goto 101
  !!-------------BndCond Drichlet Penalty Matrices----------!!
  write(9,*)"[Msg] Drichlet using Penalty method"
  tmpr1=maxval(gM1Bs1)
  tmpr2=maxval(gM1Bs4)
  !alpha=maxval(tempr(1:2))*10d5
  alpha=max(tmpr1,tmpr2)*10d5
  mass2eta=mass2
  alpha2=maxval(mass2eta)*10d5
  write(9,*)"[Info] Alpha :",alpha,"Alpha2 :",alpha2
  
  !BndCond vel eta Penalty
  do i=1,nbndpoi
    tmpi1=bndNode(i)
    k=(tmpi1-1)*ivq(0)
    if(bndNodeType(i).eq.13) then      
      if(abs(bndNormal(tmpi1,1)).gt.0.1) then
        gM1Bs1(k+ivq(tmpi1))=gM1Bs1(k+ivq(tmpi1))*alpha
      endif
      if(abs(bndNormal(tmpi1,2)).gt.0.1) then      
        gM1Bs4(k+ivq(tmpi1))=gM1Bs4(k+ivq(tmpi1))*alpha
      endif
    else if((bndNodeType(i).eq.11).or.(bndNodeType(i).eq.12).or.(bndNodeType(i).eq.14)) then
      gM1Bs1(k+ivq(tmpi1))=gM1Bs1(k+ivq(tmpi1))*alpha
      gM1Bs4(k+ivq(tmpi1))=gM1Bs4(k+ivq(tmpi1))*alpha
    endif    

    if(tmpi1.le.npoinl) then
      k=(tmpi1-1)*ivl(0)
      if((bndNodeType(i).eq.11).or.(bndNodeType(i).eq.14)) then
        mass2eta(k+ivl(tmpi1))=mass2eta(k+ivl(tmpi1)) &
          *alpha2
      endif
    endif
  enddo
  goto 102
  !!-----------End BndCond Drichlet Penalty Matrices--------!!


  !!--------------BndCond SemiDirect Matrices------------!!
  101 continue
  write(9,*)"[Msg] Drichlet using SemiDirect method"
  mass2eta=mass2
  
  !BndCond vel eta SemiDirect
  do i=1,nbndpoi
    tmpi1=bndNode(i)
    k=(tmpi1-1)*ivq(0)
    if(bndNodeType(i).eq.13) then      
      if(abs(bndNormal(tmpi1,1)).gt.0.1) then
        gM1Bs1(k+ivq(tmpi1))=1d0
        gM1Bs1(k+1:k+ivq(tmpi1)-1)=0d0
        gBs2(k+1:k+ivq(tmpi1))=0d0
      endif
      if(abs(bndNormal(tmpi1,2)).gt.0.1) then      
        gM1Bs4(k+ivq(tmpi1))=1d0
        gM1Bs4(k+1:k+ivq(tmpi1)-1)=0d0
        gBs3(k+1:k+ivq(tmpi1))=0d0
      endif
    else if((bndNodeType(i).eq.11).or.(bndNodeType(i).eq.12)&
      .or.(bndNodeType(i).eq.14))then
      gM1Bs1(k+ivq(tmpi1))=1d0
      gM1Bs1(k+1:k+ivq(tmpi1)-1)=0d0
      gBs2(k+1:k+ivq(tmpi1))=0d0
      gM1Bs4(k+ivq(tmpi1))=1d0
      gM1Bs4(k+1:k+ivq(tmpi1)-1)=0d0      
      gBs3(k+1:k+ivq(tmpi1))=0d0
    endif    

    if(tmpi1.le.npoinl) then
      k=(tmpi1-1)*ivl(0)
      if((bndNodeType(i).eq.11).or.(bndNodeType(i).eq.14)) then
        mass2eta(k+ivl(tmpi1))=1d0
        mass2eta(k+1:k+ivl(tmpi1)-1)=0d0
      endif
    endif
  enddo
  !!------------End BndCond SemiDirect Matrices----------!!

  102 continue
  ! Full Matrice A
  aFull=0
  do i=1,npoint
    k=(i-1)*ivFull(0)
    k2=(i-1)*ivq(0)
    aFull(k+1:k+ivq(i)-1)=gM1Bs1(k2+1:k2+ivq(i)-1)
    aFull(k+ivq(i):k+ivFull(i)-2)=gBs2(k2+1:k2+ivq(i)-1)
    aFull(k+ivFull(i)-1)=gBs2(k2+ivq(i))
    aFull(k+ivFull(i))=gM1Bs1(k2+ivq(i))

    k=(i+npoint-1)*ivFull(0)
    k2=(i-1)*ivq(0)
    aFull(k+1:k+ivq(i)-1)=gBs3(k2+1:k2+ivq(i)-1)
    aFull(k+ivq(i):k+ivFull(i)-2)=gM1Bs4(k2+1:k2+ivq(i)-1)
    aFull(k+ivFull(i)-1)=gBs3(k2+ivq(i))
    aFull(k+ivFull(i))=gM1Bs4(k2+ivq(i))
  enddo
  write(9,*)"[Info] afull size :",size(aFull),size(aFull)/131072d0  

  ! Deallocate A
  deallocate(gBs1,gBs2,gBs3,gBs4,poi2poi,poi2ele)

  ! Normalising matrices
  rowMaxA=0d0
  rowMaxB=0d0
  rowMaxC=0d0
  do i=1,npoinl
    k=(i-1)*ivl(0)
    rowMaxA(i)=mass2(k+ivl(i))
    rowMaxB(i)=mass2eta(k+ivl(i))
    mass2(k+1:k+ivl(i))=mass2(k+1:k+ivl(i))/rowMaxA(i)
    mass2eta(k+1:k+ivl(i))=mass2eta(k+1:k+ivl(i))/rowMaxB(i)
  enddo
  tmpi1=2*npoint
  tmpi2=2*ivq(0)
  do i=1,tmpi1
    k=(i-1)*tmpi2
    rowMaxC(i)=aFull(k+ivFull(i))
    aFull(k+1:k+tmpi2)=aFull(k+1:k+tmpi2)/rowMaxC(i)
  enddo

  !!-----------------Paralution CSR conv-----------------!!  
  ! Linear nodes
  nnzCSR1=0
  do i=1,npoinl
    nnzCSR1=nnzCSR1+ivl(i)
  enddo

  allocate(ivCSR1(npoinl+1),jvCSR1(nnzCSR1))
  allocate(mass2CSR(nnzCSR1),mass2etaCSR(nnzCSR1))

  i2=0
  ivCSR1(1)=1
  do i=1,npoinl
    ivCSR1(i+1)=ivCSR1(i)+ivl(i)
    k=(i-1)*ivl(0)
    do j=1,ivl(i)
      k2=linkl(k+j)
      i2=i2+1      
      jvCSR1(i2)=k2
      mass2CSR(i2)=mass2(k+j)
      mass2etaCSR(i2)=mass2eta(k+j)
    enddo
  enddo
  if((i2.ne.nnzCSR1).or.(ivCSR1(npoinl+1).ne.nnzCSR1+1)) then
    write(9,*)'[Err] CSR1 nnz not correct'
    stop
  endif
  deallocate(mass2eta,mass2)


  ! Quadratic nodes
  nnzCSR2=0
  do i=1,2*npoint
    nnzCSR2=nnzCSR2+ivFull(i)
  enddo

  allocate(ivCSR2(2*npoint+1),jvCSR2(nnzCSR2))
  allocate(aFullCSR(nnzCSR2))

  tmpi1=2*npoint
  i2=0
  ivCSR2(1)=1
  do i=1,tmpi1
    k=(i-1)*ivFull(0)
    ivCSR2(i+1)=ivCSR2(i)+ivFull(i)
    do j=1,ivFull(i)
      k2=linkf(k+j)
      i2=i2+1      
      jvCSR2(i2)=k2
      aFullCSR(i2)=aFull(k+j)
    enddo
  enddo  
  if((i2.ne.nnzCSR2).or.(ivCSR2(tmpi1+1).ne.nnzCSR2+1)) then
    write(9,*)'[Err] CSR2 nnz not correct'
    stop
  endif
  deallocate(aFull,ivFull,linkf)
  write(9,*)'[INFO] Solve1',npoinl,nnzCSR1
  write(9,*)'[INFO] Solve2',2*npoint,nnzCSR2

  ! Paralution init
  call paralution_init(nthrd)
  !!---------------End Paralution COO conv---------------!!  

  open(mafi(3),file="volumeData.dat")

  if(resume) then
    call resume_input(npoinl,npoint,rtime,p,q,eta,resumeFile)
  endif

  if(presOn)then
    call shipInput(probname,rTime)
    ! if(.not.resume)then
    !   call shipInitEta(npoinl,npoint,coord,eta,p,q)
    ! endif
  endif

  call outputXML(npoinl,npoint,nelem,conn,coord,probname,mafi,&
    rTime,p,q,eta,w,depth,pdt,qdt,etaMin,etaMax,porH)
  ! call output(npoinl,npoint,nelem,conn,coord,probname,mafi,&
  !   rTime,p,q,etaAbsC,w,-velAbsC,pdt,qdt,etaMin,etaMax)  

  ! Wave Probe edit
  write(9,'(a7,a,I10)')'[INF] ','Wave porbes, Num = ',probe(0)
  if(probe(0).gt.0) then    
    open(mafi(7),file="Output/waveProbes.dat")
    write(mafi(7),'(F15.6)',advance='no')rTime
    do i=1,probe(0)      
      j=probe(i)
      write(9,'(a7,2I10,2F15.6)')'[---] ',i,j,coord(j,1:2)
      write(mafi(7),'(I15,4F15.6)',advance='no')j,&
        eta(j),p(j),q(j),etadt(j)
    enddo
    write(9,*)
    write(mafi(7),*)    
  endif
  
  !!---------------------Loop Begins---------------------!!
  do tStep=1,numTSteps
    rTime=rTime+dt
    write(9,*) "Time :",tStep,":",rTime    
    call system_clock(sysClk(1))

    !Storing old vals
    wt2=wt1
    wt1=w
    etat1=eta
    pt1=p
    qt1=q
    etadtt2=etadtt1
    pdtt2=pdtt1
    qdtt2=qdtt1
    etadtt1=etadt
    pdtt1=pdt
    qdtt1=qdt
    
    
    totDe=eta+depth(1:npoinl)
    u=p(1:npoinl)/totDe
    v=q(1:npoinl)/totDe  

    if(porOn)then
      por=porV
      do i=1,npoinl      
        if(porH(i).eq.0d0)then
          por(i)=1d0
        else if(porH(i).lt.totDe(i)) then
          por(i)=((porV(i)*porH(i))+(totDe(i)-porH(i)))/totDe(i)
        endif
      enddo

      call porMatrices1(npoinl,npoint,nelem,conn,ivl,ivq,linkl,&
        linkq,invJ,depth,por,porLam,gCxF,gCyF,gAx,gAy,&
        cnst,stoneCnst)
    endif
    

    toux=cnst(9)/rhoW*(dsqrt(u**2 + v**2)*u)/(por)
    touy=cnst(9)/rhoW*(dsqrt(u**2 + v**2)*v)/(por)

    ! Dynamic Matrices
    call matrixSet2(npoinl,npoint,nelem,conn,ivl,ivq,&
      linkl,linkq,invJ,totDe,u,v,gNAdv,gGxElev,gGyElev,&
      por,porTurX,porTurY,presK,cnst,stoneCnst,presOn)    

    call bndIntegral(npoinl,npoint,nelem,nbnd,conn,mabnd,&
      invJ,depth,bndSide,bndNormal,eta,w,pdtt1,qdtt1,&
      gFp,gFq,gFn,gFw,cnst)    

    presDiv=0d0
    if(presOn) then
      call shipPosition(rTime,dt)
      call shipPress(npoint,coord,presDiv)
    endif


    !x-----------------------x!
    !!Step 1 : get w
    !w=0d0 !Allowing previous time step value as initi guess
    g4=0d0

    !g4 = D x eta  [3x3]
    do i=1,npoinl
      k=(i-1)*ivl(0)
      do j=1,ivl(i)        
        g4(i)=g4(i)+(gDMat(k+j)*etat1(linkl(k+j)))
      enddo
    enddo

    !g4 = gFw - (D x eta)
    g4=gFw-g4

    g4=g4/rowMaxA
    w=2d0*wt1-wt2    
    call paralution_fortran_solve_csr( npoinl, npoinl, nnzCSR1, &
          'GMRES' // C_NULL_CHAR, &
          'CSR' // C_NULL_CHAR, &
          'None' // C_NULL_CHAR, &
          'CSR' // C_NULL_CHAR, &
          C_LOC(ivCSR1), C_LOC(jvCSR1), &
          C_LOC(mass2CSR), C_LOC(g4), &
          errLim, 1e-15_C_DOUBLE, 1e+8_C_DOUBLE, maxIter, &
          30, 0, 1, C_LOC(w), iter, resnorm, ier )
    write(9,201)"W",ier,iter,resnorm
    if((ier.ne.1).and.(ier.ne.2)) then
      write(9,*)'[Err] Paralution error in w. ier Time ::',ier,rTime
      stop
    endif
    !x-----------------------x!

    call bndIntegral(npoinl,npoint,nelem,nbnd,conn,mabnd,&
      invJ,depth,bndSide,bndNormal,eta,w,pdtt1,qdtt1,&
      gFp,gFq,gFn,gFw,cnst)
    
    !x-----------------------x!  
    !!Step 2 : get eta
    !etadt=0d0
    g3=0d0

    !g3 = (Cx x p) + (Cy x q) [3x6] [3x6]
    do i=1,npoinl
      k=(i-1)*ivq(0)
      do j=1,ivq(i)        
        g3(i)=g3(i)+(gCxF(k+j)*pt1(linkq(k+j))) &
          +(gCyF(k+j)*qt1(linkq(k+j)))
      enddo
    enddo

    ! !g3 = Fn - ((Cx x p) + (Cy x q))
    ! g3=gFn-g3    

    ! OutletAbs eta 
    if(outAbsOn) then      
      etaAbs=0d0
      do i2=1,outAbsP(-1)
        i=outAbsP(i2)
        etaAbs(i)=etaAbsC(i)*etat1(i)
      enddo      
      ! g3=g3 - w*c3*f(x)*eta      
      do i2=1,outAbsP(-1)
        i=outAbsP(i2)
        k=(i-1)*ivl(0)
        do j=1,ivl(i)
          g3(i)=g3(i)-(mass2A(k+j)*etaAbs(linkl(k+j))) 
        enddo
      enddo      
    endif

    if(driMeth.eq.1) then
      !BndCond eta Drichlet Penalty
      do i=1,nbndpoi
        tmpi1=bndNode(i)
        if(tmpi1.le.npoinl)then
          k=(tmpi1-1)*ivl(0)
          if(bndNodeType(i).eq.11) then
            tmpr2=coord(bndNode(i),1)-waveInfo(1,1)
            tmpr3=coord(bndNode(i),2)-waveInfo(1,2)
            tmpr1=-waveA(1)*waveW(1) &
              *dcos((waveK(1)*tmpr2*waveInfo(1,4)) &
              +(waveK(1)*tmpr3*waveInfo(1,5))-(waveW(1)*rTime))
            g3(tmpi1)=rowMaxB(tmpi1)*tmpr1
          else if (bndNodeType(i).eq.14) then
            g3(tmpi1)=0d0
          endif
        endif
      enddo
    else
      !BndCond eta Drichlet SemiDirect
      do i=1,nbndpoi
        tmpi1=bndNode(i)
        if(tmpi1.le.npoinl)then      
          if(bndNodeType(i).eq.11) then
            tmpr2=coord(bndNode(i),1)-waveInfo(1,1)
            tmpr3=coord(bndNode(i),2)-waveInfo(1,2)
            tmpr1=-waveA(1)*waveW(1) &
              *dcos((waveK(1)*tmpr2*waveInfo(1,4)) &
              +(waveK(1)*tmpr3*waveInfo(1,5))-(waveW(1)*rTime))
            g3(tmpi1)=tmpr1
          else if (bndNodeType(i).eq.14) then
            g3(tmpi1)=0d0
          endif
        endif
      enddo
    endif

    g3=g3/rowMaxB
    etadt=2d0*etadtt1-etadtt2    
    call paralution_fortran_solve_csr( npoinl, npoinl, nnzCSR1, &
          'GMRES' // C_NULL_CHAR, &
          'CSR' // C_NULL_CHAR, &
          'None' // C_NULL_CHAR, &
          'CSR' // C_NULL_CHAR, &
          C_LOC(ivCSR1), C_LOC(jvCSR1), &
          C_LOC(mass2etaCSR), C_LOC(g3), &
          errLim, 1e-15_C_DOUBLE, 1e+8_C_DOUBLE, maxIter, &
          30, 0, 1, C_LOC(etadt), iter, resnorm, ier )
    write(9,201)"ETADT",ier,iter,resnorm
    if((ier.ne.1).and.(ier.ne.2)) then
      write(9,*)'[Err] Paralution error in etadt. ier Time ::',ier,rTime
      stop
    endif

    if(tStep.le.2) then
      eta=etat1+(dt*etadt)
    else
      eta=etat1+dt*(23d0*etadt-16d0*etadtt1+5d0*etadtt2)/12d0
    endif

    ! do i=1,nbndpoi
    !   tmpi1=bndNode(i)
    !   if(bndNodeType(i).eq.11) then
    !     write(9,*)"          |EtaV",rtime,eta(tmpi1)
    !     exit
    !   endif
    ! enddo    
    !x-----------------------x!

    !BndCond eta etadt Neumann  
    i2=0
    tmpr2=1d0
    tempr2(2)=1d0
    do while((tmpr2.gt.errLim).or.(tempr2(2).gt.errLim))      
      i2=i2+1
      etaRun=eta
      etadtRun=etadt
      do i=1,nbndpoi      
        if((bndNodeType(i).eq.13).or.(bndNodeType(i).eq.12))then
          tmpr1=0d0
          tempr2(1)=0d0
          tmpi1=bndNode(i)
          if(tmpi1.gt.npoinl) exit
          k=(tmpi1-1)*ivl(0)
          do j=1,(ivl(tmpi1)-1)
            tmpr1=tmpr1+(gdNdn(k+j)*etaRun(linkl(k+j)))
            tempr2(1)=tempr2(1)+(gdNdn(k+j)*etadtRun(linkl(k+j)))
          enddo            
          eta(tmpi1)=-tmpr1/gdNdn(k+ivl(tmpi1))
          etadt(tmpi1)=-tempr2(1)/gdNdn(k+ivl(tmpi1))
        endif        
      enddo    

      tmpr2=0d0
      tempr2(2)=0d0
      do i=1,nbndpoi
        if((bndNodeType(i).eq.13).or.(bndNodeType(i).eq.12)) then
          tmpr1=0d0
          tempr2(1)=0d0
          tmpi1=bndNode(i)
          if(tmpi1.gt.npoinl) exit
          k=(tmpi1-1)*ivl(0)
          do j=1,ivl(tmpi1)
            tmpr1=tmpr1+(gdNdn(k+j)*eta(linkl(k+j)))
            tempr2(1)=tempr2(1)+(gdNdn(k+j)*etadt(linkl(k+j)))                        
          enddo
          tmpr2=tmpr2+abs(tmpr1)
          tempr2(2)=tempr2(2)+abs(tempr2(1))
        endif        
      enddo
      if(i2.eq.1) then
        tmpr3=tmpr2
        tempr2(3)=tempr2(2)
        if(tmpr3.eq.0) tmpr3=1
        if(tempr2(3).eq.0) tempr2(3)=1
      endif
      tmpr2=tmpr2/tmpr3
      tempr2(2)=tempr2(2)/tmpr3
      !write(9,201)"ETA",i2,tmpr2      
      if(i2.gt.maxIter) stop
    enddo
    write(9,203)"ETA",1,i2,tmpr2,tempr2(2)    

    !Forcing Drichlet BndCond eta
    do i=1,nbndpoi
      tmpi1=bndNode(i)
      if(tmpi1.le.npoinl)then
        if(bndNodeType(i).eq.11) then
          tmpr2=coord(bndNode(i),1)-waveInfo(1,1)
          tmpr3=coord(bndNode(i),2)-waveInfo(1,2)
          tmpr1=waveA(1) &
            *dsin((waveK(1)*tmpr2*waveInfo(1,4)) &
            +(waveK(1)*tmpr3*waveInfo(1,5))-(waveW(1)*rTime))
          eta(tmpi1)=tmpr1
        endif
      endif
    enddo

    !x-----------------------x!
    !! Step 3 : get p RHS
    pdt=0d0
    g1=0d0

    !g1 = (Gx x eta) [6x3]
    do i=1,npoint
      k=(i-1)*ivl(0)
      do j=1,ivl(i)
        g1(i)=g1(i)+(gGxElev(k+j)*etat1(linkl(k+j)))      
      enddo
    enddo

    !g1 = g1 + (N x p) [6x6]
    do i=1,npoint
      k=(i-1)*ivq(0)
      do j=1,ivq(i)
        g1(i)=g1(i)+(gNAdv(k+j)*pt1(linkq(k+j))) &
              + (presK(k+j)*presDiv(linkq(k+j),1))
      enddo
    enddo

    !g1 = g1 + (Ax x w) + (T x toux) [6x3]
    do i=1,npoint
      k=(i-1)*ivl(0)
      do j=1,ivl(i)
        g1(i)=g1(i)+(gAx(k+j)*w(linkl(k+j)))
        g1(i)=g1(i)+(touMass(k+j)*toux(linkl(k+j)))
        !g1(i)=g1(i)+(touMass(k+j)*presDiv(linkl(k+j),1))
      enddo
    enddo

    !g1 = Fp - ((N x p) + (Gx x eta)) - (Ax x w)    
    g1=gFp-g1    

    ! Porosity
    !g1=g1 - (porLam*p) - (porTurX*|p|) [6x6] [6x6]
    if(porOn) then
      do i=1,npoint
        k=(i-1)*ivq(0)
        do j=1,ivq(i)
          ! g1(i)=g1(i)-(porLam(k+j)*pt1(linkq(k+j))) &
          !   -(porTurX(k+j)*abs(pt1(linkq(k+j))))
          g1(i)=g1(i)-(porLam(k+j)*pt1(linkq(k+j))) &
            -(porTurX(k+j)*dsqrt(pt1(linkq(k+j))**2 &
            + qt1(linkq(k+j))**2))
        enddo
      enddo
    endif
    !x-----------------------x!

    !x-----------------------x!
    !! Step 4 : get q RHS
    qdt=0d0
    g2=0d0

    !g2 = (Gy x eta) [6x3]
    do i=1,npoint
      k=(i-1)*ivl(0)
      do j=1,ivl(i)
        g2(i)=g2(i)+(gGyElev(k+j)*etat1(linkl(k+j)))
        !g2(i)=g2(i)+(-gAy(k+j)/cnst(2)*etat1(linkl(k+j)))
      enddo
    enddo

    !g2 = g2 + (N x q) [6x6]
    do i=1,npoint
      k=(i-1)*ivq(0)
      do j=1,ivq(i)
        g2(i)=g2(i)+(gNAdv(k+j)*qt1(linkq(k+j))) &
              + (presK(k+j)*presDiv(linkq(k+j),2))
      enddo
    enddo

    !g2 = g2 + (Ay x w) + (T x touy) [6x3]
    do i=1,npoint
      k=(i-1)*ivl(0)
      do j=1,ivl(i)
        g2(i)=g2(i)+(gAy(k+j)*w(linkl(k+j)))
        g2(i)=g2(i)+(touMass(k+j)*touy(linkl(k+j)))
        !g2(i)=g2(i)+(touMass(k+j)*presDiv(linkl(k+j),2))
      enddo
    enddo

    !g2 = Fq - ((N x q) + (Gy x eta)) - (Ay x w)    
    g2=gFq-g2    

    ! Porosity
    !g2=g2 - (porLam*q) - (porTurY*q) [6x6] [6x6]
    if(porOn) then
      do i=1,npoint
        k=(i-1)*ivq(0)
        do j=1,ivq(i)
          ! g2(i)=g2(i)-(porLam(k+j)*qt1(linkq(k+j))) &
          !   -(porTurY(k+j)*abs(qt1(linkq(k+j))))
          g2(i)=g2(i)-(porLam(k+j)*qt1(linkq(k+j))) &
            -(porTurY(k+j)*dsqrt(pt1(linkq(k+j))**2 &
            + qt1(linkq(k+j))**2))
        enddo
      enddo
    endif
    !x-----------------------x!
    

    ! OutletAbs vel
    if(outAbsOn) then      
      velAbs=0d0
      do i2=1,outAbsP(0)
        i=outAbsP(i2)
        velAbs(i,1)=velAbsC(i)*pt1(i)
        velAbs(i,2)=velAbsC(i)*qt1(i)
      enddo      
      ! g1=g1 - w*c1*f(x)*u
      do i2=1,outAbsP(0)
        i=outAbsP(i2)
        k=(i-1)*ivq(0)
        do j=1,ivq(i)
          g1(i)=g1(i)-(mass1(k+j)*velAbs(linkq(k+j),1))
          g2(i)=g2(i)-(mass1(k+j)*velAbs(linkq(k+j),2))
        enddo
      enddo      
    endif

    ! if(outAbsOn) then
    !   ! OutletAbs vel
    !   velAbs=0d0
    !   do i=1,npoint
    !     if(coord(i,1).gt.outAbs(1)) then
    !       tmpr1=coord(i,1)-outAbs(1)
    !       tmpr2=outAbs(4) &
    !         *(exp((tmpr1/outAbs(3))**2)-1d0)/outAbs(5)
    !       velAbs(i,1)=tmpr2*pt1(i)
    !       velAbs(i,2)=tmpr2*qt1(i)
    !     endif
    !   enddo
    !   ! g1=g1 - w*c1*f(x)*u
    !   do i=1,npoint
    !     if(coord(i,1).gt.outAbs(1)) then
    !       k=(i-1)*ivq(0)
    !       do j=1,ivq(i)
    !         g1(i)=g1(i)-(mass1(k+j)*velAbs(linkq(k+j),1))
    !         g2(i)=g2(i)-(mass1(k+j)*velAbs(linkq(k+j),2))
    !       enddo
    !     endif
    !   enddo    
    ! endif
 
    !x-----------------------x!
    !! Step 5 : get p and q

    if(driMeth.eq.1)then
      !BndCond vel Drichlet Penalty
      do i=1,nbndpoi
        tmpi1=bndNode(i)
        if(bndNodeType(i).eq.13)then
          if(abs(bndNormal(tmpi1,1)).gt.0.1) then
            g1(tmpi1)=0
          endif
          if(abs(bndNormal(tmpi1,2)).gt.0.1) then
            g2(tmpi1)=0
          endif        
        endif

        if((bndNodeType(i).eq.12).or.(bndNodeType(i).eq.14))then        
          g1(tmpi1)=0        
          g2(tmpi1)=0
        endif

        if(bndNodeType(i).eq.11)then
          k=(tmpi1-1)*ivq(0)
          tmpr3=coord(bndNode(i),1)-waveInfo(1,1)
          tmpr4=coord(bndNode(i),2)-waveInfo(1,2)
          tmpr1=-waveA(1)*waveW(1)*waveW(1)/waveK(1)&
            *waveInfo(1,4)*dcos((waveK(1)*tmpr3*waveInfo(1,4)) &
            +(waveK(1)*tmpr4*waveInfo(1,5))-(waveW(1)*rTime))
          tmpr2=-waveA(1)*waveW(1)*waveW(1)/waveK(1)&
            *waveInfo(1,5)*dcos((waveK(1)*tmpr3*waveInfo(1,4)) &
            +(waveK(1)*tmpr4*waveInfo(1,5))-(waveW(1)*rTime))
          g1(tmpi1)=rowMaxC(tmpi1)*tmpr1                
          g2(tmpi1)=rowMaxC(npoint+tmpi1)*tmpr2
        endif
      enddo                    
    else
      !BndCond vel Drichlet SemiDirect
      do i=1,nbndpoi
        tmpi1=bndNode(i)
        if(bndNodeType(i).eq.13)then
          if(abs(bndNormal(tmpi1,1)).gt.0.1) then
            g1(tmpi1)=0
          endif
          if(abs(bndNormal(tmpi1,2)).gt.0.1) then
            g2(tmpi1)=0
          endif        
        endif

        if((bndNodeType(i).eq.12).or.(bndNodeType(i).eq.14))then        
          g1(tmpi1)=0        
          g2(tmpi1)=0
        endif

        if(bndNodeType(i).eq.11)then        
          tmpr3=coord(bndNode(i),1)-waveInfo(1,1)
          tmpr4=coord(bndNode(i),2)-waveInfo(1,2)
          tmpr1=-waveA(1)*waveW(1)*waveW(1)/waveK(1)&
            *waveInfo(1,4)*dcos((waveK(1)*tmpr3*waveInfo(1,4)) &
            +(waveK(1)*tmpr4*waveInfo(1,5))-(waveW(1)*rTime))
          tmpr2=-waveA(1)*waveW(1)*waveW(1)/waveK(1)&
            *waveInfo(1,5)*dcos((waveK(1)*tmpr3*waveInfo(1,4)) &
            +(waveK(1)*tmpr4*waveInfo(1,5))-(waveW(1)*rTime))
          g1(tmpi1)=tmpr1                
          g2(tmpi1)=tmpr2
        endif
      enddo                    
    endif

    xFull(1:npoint)=2d0*pdtt1-pdtt2
    xFull(npoint+1:2*npoint)=2d0*qdtt1-qdtt2
    bFull(1:npoint)=g1
    bFull(npoint+1:2*npoint)=g2
    
    bFull=bFull/rowMaxC
    call system_clock(sysClk(2))
    tmpi1=2*npoint
    call paralution_fortran_solve_csr( tmpi1, tmpi1, nnzCSR2, &
          'GMRES' // C_NULL_CHAR, &
          'CSR' // C_NULL_CHAR, &
          'None' // C_NULL_CHAR, &
          'CSR' // C_NULL_CHAR, &
          C_LOC(ivCSR2), C_LOC(jvCSR2), &
          C_LOC(aFullCSR), C_LOC(bFull), &
          errLim, 1e-15_C_DOUBLE, 1e+8_C_DOUBLE, maxIter, &
          30, 0, 1, C_LOC(xFull), iter, resnorm, ier )
    write(9,201)"P Q",ier,iter,resnorm
    if((ier.ne.1).and.(ier.ne.2)) then
      write(9,*)'[Err] Paralution error in P Q. ier Time ::',ier,rTime
      stop
    endif
    call system_clock(sysClk(3))

    
    pdt=xFull(1:npoint)
    qdt=xFull(npoint+1:2*npoint)    


    ! Time stepping
    if(tStep.le.2) then
      p=pt1+(pdt*dt)
      q=qt1+(qdt*dt)
    else 
      p=pt1+dt*(23d0*pdt-16d0*pdtt1+5d0*pdtt2)/12d0
      q=qt1+dt*(23d0*qdt-16d0*qdtt1+5d0*qdtt2)/12d0
    endif

    !Forcing Drichlet BndCond vel
    do i=1,nbndpoi
      tmpi1=bndNode(i)      
      if(bndNodeType(i).eq.11) then
        tmpr3=coord(bndNode(i),1)-waveInfo(1,1)
        tmpr4=coord(bndNode(i),2)-waveInfo(1,2)
        tmpr1=waveA(1)*waveW(1)/waveK(1)&
          *waveInfo(1,4)*dsin((waveK(1)*tmpr3*waveInfo(1,4)) &
          +(waveK(1)*tmpr4*waveInfo(1,5))-(waveW(1)*rTime))
        tmpr2=waveA(1)*waveW(1)/waveK(1)&
          *waveInfo(1,5)*dsin((waveK(1)*tmpr3*waveInfo(1,4)) &
          +(waveK(1)*tmpr4*waveInfo(1,5))-(waveW(1)*rTime))
        p(tmpi1)=tmpr1                
        q(tmpi1)=tmpr2
      endif      
    enddo

    tmpr1=min(mod(rTime,5d0),abs(5d0-mod(rTime,5d0)))
    if(tmpr1.lt.dt/1e1) then
      write(9,'("      |",a5,a)')'[INF]',&
        ' Resetting waveHeight'
      etaMin=0d0
      etaMax=0d0
    endif    
    do i=1,npoinl
      if(eta(i).lt.etaMin(i)) etaMin(i)=eta(i)
      if(eta(i).gt.etaMax(i)) etaMax(i)=eta(i)
    enddo
    !x-----------------------x!

    ! Wave Probe edit
    if(probe(0).gt.0) then
      write(mafi(7),'(F15.6)',advance='no')rTime
      do i=1,probe(0)
        j=probe(i)
        write(mafi(7),'(I15,4F15.6)',advance='no')j,&
          eta(j),p(j),q(j),etadt(j)
      enddo
      write(mafi(7),*)
    endif

    ! Volume Check
    tmpr1=0d0
    totDe=eta+depth(1:npoinl)
    do i=1,nelem
      nl(1:3)=conn(i,1:3)
      tmpr2=0d0
      do j=1,3
        tmpr2=tmpr2+totDe(nl(j))
      enddo
      tmpr2=tmpr2/6d0*invJ(i,5)
      tmpr1=tmpr1+tmpr2
    enddo

    call system_clock(sysClk(4))
    tmpr1=(sysClk(4)-sysClk(1))/sysRate
    tmpr2=(sysClk(3)-sysClk(2))/sysRate
    write(9,202)"[SPD]",tmpr1,tmpr2,tmpr2/tmpr1*100d0
      

    if(mod(tStep,fileOut).eq.0) then
      write(mafi(3),*)tStep*dt,tmpr1
      call outputXML(npoinl,npoint,nelem,conn,coord,probname,&
        mafi,rTime,p,q,eta,w,depth,pdt,qdt,etaMin,etaMax,&
        porH)
    endif

    if(mod(tStep,resumeOut).eq.0) then      
      call resume_output(npoinl,npoint,probname,mafi,&
        rTime,p,q,eta)
    endif

  enddo
  !! Loop Ends
  close(mafi(3))
  if(probe(0).gt.0) close(mafi(7))

  201 format('      |',a5,i10,i10,e15.4)
  202 format('      |',a5,3f15.4)
  203 format('      |',a5,i10,i10,2e15.4)

  ! Paralution backend stop
  call paralution_stop
  
  !!-----------------------End Loop----------------------!!
!!-------------------------End Loop---------------------------!!

  write(9,*)"boussinesqQuad End"
end program boussinesqQuad
