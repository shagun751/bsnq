!!---------------------- Version 1.x.x -----------------------!!
!!  -> Quadratic + Linear
!!  -> Quadratic Jacobian
!!  -> Superparametric - eta, w
!!  -> Isoparametric - P Q
!!  -> Boundary - SemiDirect + Penalty + Gauss Seidel
!!    -> 11 - Inlet - No Absorbing
!!    -> 12 - NoSlip Wall 
!!    -> 13 - Slip Wall - Rectangular wall only
!!    -> 14 - Outlet  - Not coded
!!    -> 15 - Sponge - BOUSS2D approach generalised input
!!  -> Porosity - Emergent structure only
!!    -> Generalised input
!!  -> Bottom Shear
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
!!---------------------- Version 1.x.x -----------------------!!
!!  Mar 16 2018
!!  1000
!!  continued from bsnq_v8.36

include 'subroutines/mods.f90'
!include 'subroutines/bndIntegral.f90'
!include 'subroutines/boundaryDifferential.f90'
include 'subroutines/geometry.f90'
! include 'subroutines/gqMatrixSet1.f90'
! include 'subroutines/matrixSet1.f90'
! include 'subroutines/matrixSet2.f90'
include 'subroutines/mergeSort.f90'
include 'subroutines/nodeConnAll.f90'
include 'subroutines/outputN.f90'
! include 'subroutines/porMatrices.f90'
! include 'subroutines/resumeFile.f90'
! include 'subroutines/shipFncs.f90'
! include 'subroutines/waveCalculator_v2.f90'


program boussinesqQuad

use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE
use, intrinsic :: ISO_C_BINDING, only : C_CHAR, C_NULL_CHAR, C_LOC
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

!!--------------------------Declarations----------------------!!  
  integer(kind=C_K1)::npl,npq,npt,nele,nedg,nbndtyp,nbnd
  integer(kind=C_K1),allocatable::conn(:,:),mabnd(:,:),tempia(:,:)
  integer(kind=C_K1),allocatable::poi2poi(:,:),poi2ele(:,:)
  integer(kind=C_K1),allocatable::npoisur(:,:),bndNode(:,:)
  integer(kind=C_K1),allocatable::tmpia(:,:),probe(:)
  integer(kind=C_K1),allocatable::outAbsN(:)
  integer(kind=C_K1)::numTSteps,tStep,ier,maxIter,iter,gsIter
  integer(kind=C_K1)::nthrd,driMeth,fileOut,resumeOut
  integer(kind=C_K1)::nnzCSR1,nnzCSR2

  integer(kind=C_K1),allocatable::ivl(:),linkl(:),ivq(:),linkq(:)
  integer(kind=C_K1),allocatable::ivFull(:),linkf(:)

  real(kind=C_K2)::dt,errLim,rTime
  real(kind=C_K2)::wavT,wavH,wavCx,wavCy,wavTh,wavSn,wavCs
  real(kind=C_K2),allocatable::coor(:,:),dep(:),tmpra(:,:)
  real(kind=C_K2),allocatable::bndEdg(:,:),nodNor(:,:) 
  real(kind=C_K2),allocatable::w(:),eta(:),etat1(:)
  real(kind=C_K2),allocatable::tDe(:),tDet1(:)  
  real(kind=C_K2),allocatable::p(:),pt1(:),q(:),qt1(:)
  real(kind=C_K2),allocatable::edtt1(:),pdtt1(:),qdtt1(:)  
  real(kind=C_K2),allocatable::por(:),porH(:),porV(:)
  real(kind=C_K2),allocatable::outAbsV(:,:)

  type(jacbType),allocatable::jacb(:)

  character(len=100)::probname,text,resumeFile

  logical::ex,resume,presOn,outAbsOn,porOn  
  logical,allocatable::templa(:)       
!!------------------------End Declarations--------------------!!

!!--------------------------All Constants---------------------!!

  call getarg(1,probname)  
  if(len_trim(probname).lt.1) then
    write(*,'(A)',advance='no')"Enter Problem Name:"
    read(*,*)probname
  endif
  write(*,*)"[MSG] Problem Name: "//probname(1:len_trim(probname))
!!------------------------End All Constants-------------------!!  

!!-------------------------Mesh File Reading------------------!!
  !Opening mesh file  
  text=probname(1:len_trim(probname))//'.plt'
  inquire(file=text(1:len_trim(text)),exist=ex)
  if(ex) then
    open(mafi(1),file=text(1:len_trim(text)))    
  else
    write(*,*)"[ERR] Missing mesh file"
    stop
  endif

  text=probname(1:len_trim(probname))//'.rout'
  open(9,file=text(1:len_trim(text)))

  !!------------------------Mesh Type 0------------------!!
  read(mafi(1),*,end=11,err=11)text
  read(mafi(1),*,end=11,err=11)nele,npl,nedg
  read(mafi(1),*,end=11,err=11)text
  read(mafi(1),*,end=11,err=11)nbnd,nbndtyp
  write(9,'(" [MSG] ",3a15)')"Elements","Linear Nodes","Edges"
  write(9,'(" [---] ",3i15)')nele,npl,nedg
  write(9,'(" [---] ",2a15)')"Bnd","BndTypes"
  write(9,'(" [---] ",2i15)')nbnd,nbndtyp

  !Assumption regarding number of quad Nodes
  npq=nedg
  npt=npl+npq

  !Nodes Read
  allocate(coor(npt,2))
  coor=-999  
  read(mafi(1),*,end=11,err=11)text
  do i=1,npl
    read(mafi(1),*,end=11,err=11)coor(i,1),coor(i,2)
  enddo
  write(9,*)"[MSG] Nodes read done"
  !Debug comments
  ! do i=1,npoinl
  !   write(9,*)i,":",(coord(i,j),j=1,2),depth(i)
  ! enddo

  !Elements Read
  allocate(conn(nele,6))
  conn=0
  read(mafi(1),*,end=11,err=11)text
  do i=1,nele
    read(mafi(1),*,end=11,err=11)conn(i,3),conn(i,1),conn(i,2)    
  enddo
  write(9,*)"[MSG] Elements read done"
  !Debug comments
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
  write(9,*)"[MSG] Boundaries read done" 
  ! ! Debug comments 
  ! do i=1,nbnd
  !   write(9,*)i,":",(mabnd(i,j),j=1,6)
  ! enddo

  !Depth Read
  allocate(dep(npt))
  dep=-999
  read(mafi(1),*,end=11,err=11)text
  do i=1,npl
    read(mafi(1),*,end=11,err=11)dep(i)
  enddo
  write(9,*)"[MSG] Depth read done"
  !Debug comments
  ! do i=1,npoinl
  !   write(9,*)i,depth(i)
  ! enddo

  goto 12
  !!----------------------End Mesh Type 0----------------!!

  11 write(9,*) "[ERR] Check mesh file format"
  stop
  13 write(9,*) "[ERR] hex2dec error"
  stop
  14 write(9,*) "[ERR] Number of boundaries mismatch"
  stop
  12 close(mafi(1))
  write(9,*) "[MSG] Mesh file read successful"
  write(9,*)
!!---------------------End Mesh File Reading------------------!!

!!--------------------Generate Middle Points------------------!!
  allocate(poi2poi(npt,maxNePoi),poi2ele(npt,maxNeEle))
  allocate(npoisur(npt,3))
  poi2poi=0
  poi2ele=0
  npoisur=0

  call middleNode(npl,npq,npt,nele,nedg,maxNePoi,&
    coor,dep,conn,poi2poi,npoisur)
  if(npq.ne.nedg) then
    write(9,*) "[ERR] Initial npoint assumption insufficiant"
    stop
  endif
  poi2poi=0
  npoisur=0

  !!Debug comments
  ! do i=1,npt
  !   write(9,'(i10," : ",3f15.6)')i,coor(i,:),dep(i)
  ! enddo
  ! do i=1,nele
  !   write(9,'(i10," : ",6i10)')i,(conn(i,j),j=1,6)
  ! enddo

  write(9,*)"[MSG] Middle points generation done"
  write(9,'(" [MSG] ",3a15)')"Lin Nod","Quad Nod","Tot Nod"
  write(9,'(" [---] ",3i15)')npl,npq,npt
  write(9,'(" [---] ",3a15)')"Elements","Bnd","BndTypes"
  write(9,'(" [---] ",3i15)')nele,nbnd,nbndtyp

  !! Boundary sides middle nodes update
  !! middle nodes stored in mabnd(i,5)
  call bndMiddleNodes(npl,npt,nele,nbnd,conn,mabnd)
  write(9,*) "[MSG] Boundary sides middle point done"  
  !! Debug comments
  ! do i=1,nbnd
  !   write(9,'(i10,":",6i10)')i,mabnd(i,:)
  ! enddo

  !! Storing boundary nodes
  tmpi5=0
  allocate(tmpia(npt,2))
  tmpia=0
  do i=1,nbnd
    nq(1)=mabnd(i,1)
    nq(2)=mabnd(i,2)
    nq(3)=mabnd(i,5)
    tmpi4=mabnd(i,4)
    do j=1,3
      do k=1,tmpi5
        if(tmpia(k,1).eq.nq(j)) goto 31
      enddo
      tmpi5=tmpi5+1
      tmpia(tmpi5,1)=nq(j)
      tmpia(tmpi5,2)=tmpi4
      k=tmpi5
      31 continue
      !Preferance for bndtype 14 over 12 or 13
      if((tmpia(k,2).ne.14).and.(tmpi4.eq.14)) then
        tmpia(k,2)=14
      end if
      !Preferance for bndtype 11 - higher
      if((tmpia(k,2).ne.11).and.(tmpi4.eq.11)) then
        tmpia(k,2)=11
      end if
      ! if((tempia(k,2).eq.12).and.(tmpi4.eq.13)) then
      !   tempia(k,2)=13
      ! end if
    enddo
  enddo
  !! bndNode Information
  !! bndNode(0,0)   = Total Num of bnd Nodes
  !! bndNode(-1,0)  = Index of last 11 node
  !! bndNode(-2,0)  = Index of last 12 node
  !! bndNode(-3,0)  = Index of last 13 node
  !! bndNode(-4,0)  = Index of last 14 node
  allocate(bndNode(-4:tmpi5,2))
  bndNode=0
  bndNode(0,1)=tmpi5
  k=0
  do i=11,14
    do j=1,tmpi5
      if(tmpia(j,2).ne.i)cycle
      k=k+1 
      bndNode(k,:)=tmpia(j,:)
    enddo    
    bndNode(10-i,1)=k
  enddo
  deallocate(tmpia)
  if(bndNode(0,1).ne.bndNode(-4,1))then
    write(9,*)"[ERR] Check boundary node types more than typ 14"
  endif
  ! !! Debug comments    
  ! do i=-4,bndNode(0,1)    
  !   write(9,'(i10," : ",2i10)')i,bndNode(i,:)
  ! enddo
  write(9,*)"[MSG] Storing bndNodes done"
  write(9,*)
!!------------------End Generate Middle Points----------------!!

!!-----------------------Node Connectivity--------------------!!
  poi2poi=0
  poi2ele=0
  npoisur=0
  call nodeConnVSR(npt,nele,maxNePoi,maxNeEle,conn,&
    poi2poi,poi2ele,npoisur)
  do i=1,npt
    call mergeSort(maxNePoi,npoisur(i,1),poi2poi(i,:))
  enddo

  !! Finding number of linear element nbhs
  do i=1,npt
    do j=1,npoisur(i,1)
      if(poi2poi(i,j).gt.npl) exit
      npoisur(i,3)=npoisur(i,3)+1
    enddo
  enddo

  ! !!Debug Comments
  ! do i=1,npt
  !   write(9,*) i,":",(poi2poi(i,j),j=1,maxNePoi)
  ! enddo
  ! do i=1,npt
  !   write(9,*) i,":",(poi2ele(i,j),j=1,maxNeEle)
  ! enddo
  ! do i=1,npt
  !   write(9,*) i,":",(npoisur(i,j),j=1,3)
  ! enddo
  write(9,*)"[MSG] Node Connectivity Done"
!!---------------------End Node Connectivity------------------!!

!!---------------------------VSR Matrices---------------------!!
  allocate(ivl(0:npt),ivq(0:npt))  
  ivl(0)=maxval(npoisur(:,3))+1
  ivq(0)=maxval(npoisur(:,1))+1
  allocate(linkl(ivl(0)*npt),linkq(ivq(0)*npt))
  linkl=0
  linkq=0

  !! IV matrix linear and quadratic
  do i=1,npl
    ivl(i)=npoisur(i,3)+1
  enddo
  do i=npl+1,npt
    ivl(i)=npoisur(i,3)
  enddo
  ! write(9,*) ivl

  do i=1,npt
    ivq(i)=npoisur(i,1)+1
  enddo
  ! write(9,*) ivq

  !! Link matrix linear
  do i=1,npl
    k=(i-1)*ivl(0)
    do j=1,npoisur(i,3)
      linkl(k+j)=poi2poi(i,j)
    enddo
    linkl(k+ivl(i))=i
  enddo
  do i=npl+1,npt
    k=(i-1)*ivl(0)
    do j=1,npoisur(i,3)
      linkl(k+j)=poi2poi(i,j)
    enddo    
  enddo
  ! do i=1,npt
  !   k=(i-1)*ivl(0)
  !   write(9,*) i,":",linkl(k+1:k+ivl(i))
  ! enddo

  !! Link matrix quadratic
  do  i=1,npt
    k=(i-1)*ivq(0)
    do j=1,npoisur(i,1)
      linkq(k+j)=poi2poi(i,j)
    enddo
    linkq(k+ivq(i))=i
  enddo
  ! write(9,*)linkq

  !! Full Matrices ivFull and linkf
  allocate(ivFull(0:2*npt))
  ivFull(0)=2*ivq(0)
  allocate(linkf(ivFull(0)*2*npt))
  do i=1,npt
    ivFull(i)=2*ivq(i)
    ivFull(npt+i)=2*ivq(i)
  enddo

  do i=1,npt
    k=(i-1)*ivFull(0)
    linkf(k+1:k+npoisur(i,1))=poi2poi(i,1:npoisur(i,1))
    linkf(k+npoisur(i,1)+1:k+(2*npoisur(i,1)))=npt &
      + poi2poi(i,1:npoisur(i,1))
    linkf(k+ivFull(i)-1)=npt+i
    linkf(k+ivFull(i))=i

    k=(i+npt-1)*ivFull(0)
    linkf(k+1:k+npoisur(i,1))=poi2poi(i,1:npoisur(i,1))
    linkf(k+npoisur(i,1)+1:k+(2*npoisur(i,1)))=npt &
      + poi2poi(i,1:npoisur(i,1))
    linkf(k+ivFull(i)-1)=i
    linkf(k+ivFull(i))=npt+i
  enddo
  write(9,*)"[MSG] VSR storage matrices done"
  write(9,*)
!!-------------------------End VSR Matrices-------------------!!

!!-------------------------Jacobian and Normals---------------!!
  allocate(jacb(nele),bndEdg(nbnd,3),nodNor(npt,2))
  call initGaussPoi
  call initShapeFnc
  call jacbInvQuad(npt,nele,conn,coor,jacb)
  !! Priniting area using quad jacb
  ! do iel=1,nele    
  !   tmpr2=0d0
  !   do i=1,nGP
  !     tmpr2=tmpr2+gpW(i)*jacb(iel)%D(i)
  !   enddo
  !   !if(abs(tmpr2-tmpr1).lt.1e-10)cycle
  !   write(*,'(i10,f15.8)')iel,tmpr2
  ! enddo  
  ! stop

  !! Boundary side normals
  call bndSideInfo(npl,npt,nele,nbnd,coor,mabnd,bndEdg)
  ! do i=1,nbnd
  !   write(9,*)i,":",bndEdg(i,:),mabnd(i,:)
  ! enddo

  !! Boundary Node Normal
  nodNor=0d0
  call bndNodeNormal(npl,npt,nele,nbnd,coor,mabnd,bndEdg,nodNor)
  ! do i=1,npt
  !   write(9,'(i10," : ",2f15.6)')i,nodNor(i,:)
  ! enddo
!!-----------------------End Jacobian and Normals-------------!!

!!------------------------Allocations-------------------------!!
  allocate(w(npl),eta(npl),etat1(npl))
  allocate(tDe(npt),tDet1(npt))
  allocate(por(npt),porH(npt),porV(npt))
  allocate(p(npt),pt1(npt),q(npt),qt1(npt))
  allocate(edtt1(npl),pdtt1(npt),qdtt1(npt))
!!----------------------End Allocations-----------------------!!

!!--------------------------Settings--------------------------!!
  !!-------------------------Input-----------------------!!
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

  !! Inlet Wave Characteristics
  read(mafi(5),*,end=81,err=81)text
  read(mafi(5),*,end=81,err=81)wavT,wavH
  read(mafi(5),*,end=81,err=81)text
  read(mafi(5),*,end=81,err=81)wavCx,wavCy,wavTh
  wavTh=wavTh*pi/180d0
  wavSn=dsin(wavTh)
  wavCs=dcos(wavTh)

  !! OutletAbs Coefficient
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
    outAbsV(:,3)=30d0/wavT
    outAbsV(:,4)=exp(1d0)-1d0
  endif

  !! Porosity
  por=1d0
  porV=1d0
  porH=0d0
  read(mafi(5),*,end=81,err=81)text
  read(mafi(5),*,end=81,err=81)text
  read(mafi(5),*,end=81,err=81)porOn    
  if(porOn)then      
    text=probname(1:len_trim(probname))//'_por.plt'
    inquire(file=text(1:len_trim(text)),exist=ex)
    if(ex) then
      open(mafi(6),file=text(1:len_trim(text)))        
    else
      write(9,*)"[ERR] Missing porosity value file"
      stop
    endif

    read(mafi(6),*,end=91,err=91)text
    do i=1,npl
      read(mafi(6),*,end=91,err=91)porV(i),porH(i)
    enddo

    goto 92
    91 write(9,*) "[ERR] Check porosity file format"
    stop
    92 close(mafi(6))
    write(9,*)'[MSG] Porosity file read'  
  endif

  goto 82
  81 write(9,*) "[ERR] Check input file format"
  stop
  82 close(mafi(5))
  !!-----------------------End Input---------------------!!

!!------------------------End Settings------------------------!!
  rTime=0d0
  w=0d0
  eta=0d0;  etat1=0d0;
  p=0d0;  pt1=0d0;
  q=0d0;  qt1=0d0;
  edtt1=0d0;  pdtt1=0d0;  qdtt1=0d0;


!!-------------------------Initialise-------------------------!!

  call outputXML(npl,npt,nele,conn,coor,probname,mafi,&
    rTime,p,q,eta,w,dep)
!!-------------------------Initialise-------------------------!!

  write(9,*)"boussinesqQuad End"
end program boussinesqQuad
