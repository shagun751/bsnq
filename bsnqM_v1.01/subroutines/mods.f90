!!----------------------------basicVars----------------------------!!
module bsnqGlobVars 
use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE
use, intrinsic :: ISO_C_BINDING, only : C_CHAR, C_NULL_CHAR, C_LOC 
implicit none
  
  integer,parameter::C_K1=C_INT,C_K2=C_DOUBLE,C_KSTR=256,C_LG=1  

  !!-----------------------Constants-----------------------!!
  !! Gravity
  real(kind=C_K2),parameter::grav=9.81d0          
  !! PI
  real(kind=C_K2),parameter::pi=atan(1d0)*4d0       
  real(kind=C_K2),parameter::deg2rad=pi/180d0
  real(kind=C_K2),parameter::rad2deg=180d0/pi
  !! Water density
  real(kind=C_K2),parameter::rhoW=1000d0  
  !! Bsnq constant
  real(kind=C_K2),parameter::BsqC=1d0/15d0
  !!---------------------End Constants---------------------!!


  !!  mafi defintions
  !!  mafi(1)     Mesh File
  !!  mafi(2)     Paraview output
  !!  mafi(3)     Volume output
  !!  mafi(4)     <Unknown>
  !!  mafi(5)     Input file
  !!  mafi(6)     Porosity file
  !!  mafi(7)     Wave probes files


end module bsnqGlobVars
!!--------------------------End basicVars--------------------------!!



!!-------------------------airyWaveModule--------------------------!!
module airyWaveModule
use bsnqGlobVars
implicit none

  type, public :: airyType    
    real(kind=C_K2)::T,d,H,L,w
    real(kind=C_K2)::x0,y0
    real(kind=C_K2)::k,kx,ky
    real(kind=C_K2)::thDeg,thRad,csth,snth      
  contains
    procedure ::  getEta
    procedure ::  getPQ
  end type airyType

  interface airyType
     procedure :: waveLenCalc
  end interface airyType

contains

!!-------------------------waveLenCalc-------------------------!!
  type(airyType) function waveLenCalc(inT,inD,inH,inX0,inY0,inThDeg)
  implicit none

    !! WaveLen using dispersion relation from Airy wave-theory

    integer(kind=C_K1)::iterMax,i
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
    waveLenCalc%kx=waveLenCalc%k*waveLenCalc%csth
    waveLenCalc%ky=waveLenCalc%k*waveLenCalc%snth

  end function waveLenCalc
!!-----------------------End waveLenCalc-----------------------!!



!!---------------------------getEta----------------------------!!
  subroutine getEta(b,rTime,x,y,eta)
  implicit none

    class(airyType),intent(in)::b

    !integer(kind=C_K1)::

    real(kind=C_K2),intent(in)::rTime,x,y
    real(kind=C_K2),intent(out)::eta
    real(kind=C_K2)::dx,dy

    dx=(x-b%x0)
    dy=(y-b%y0)
    eta=b%H/2d0 * dsin(b%kx*dx + b%ky*dy - b%w*rTime)

  end subroutine getEta
!!-------------------------End getEta--------------------------!!



!!----------------------------getPQ----------------------------!!
  subroutine getPQ(b,rTime,x,y,p,q)
  implicit none

    class(airyType),intent(in)::b

    !integer(kind=C_K1)::

    real(kind=C_K2),intent(in)::rTime,x,y
    real(kind=C_K2),intent(out)::p,q
    real(kind=C_K2)::dx,dy

    dx=(x-b%x0)
    dy=(y-b%y0)
    p=b%H/2d0 * b%w / b%k * b%csth &
      * dsin(b%kx*dx + b%ky*dy - b%w*rTime)
    q=b%H/2d0 * b%w / b%k * b%snth &
      * dsin(b%kx*dx + b%ky*dy - b%w*rTime)

  end subroutine getPQ
!!--------------------------End getPQ--------------------------!!



!!-------------------------getEtadxdy--------------------------!!
  subroutine getEtadxdy(b,rTime,x,y,detadx,detady)
  implicit none

    class(airyType),intent(in)::b

    !integer(kind=C_K1)::

    real(kind=C_K2),intent(in)::rTime,x,y
    real(kind=C_K2),intent(out)::detadx,detady
    real(kind=C_K2)::dx,dy

    dx=(x-b%x0)
    dy=(y-b%y0)
    detadx=b%H/2d0 * b%kx * dcos(b%kx*dx + b%ky*dy - b%w*rTime)
    detady=b%H/2d0 * b%ky * dcos(b%kx*dx + b%ky*dy - b%w*rTime)

  end subroutine getEtadxdy
!!-----------------------End getEtadxdy------------------------!!

end module airyWaveModule
!!-----------------------End airyWaveModule------------------------!!



!!--------------------------outAbsModule---------------------------!!
module outAbsModule
use bsnqGlobVars
implicit none

  type, public :: absTyp
    integer(kind=C_K1)::N,typ
    real(kind=C_K2)::x0,l,c1,wT
  contains
    procedure ::  calcAbsC
  end type absTyp

  interface absTyp
    procedure ::  initAbsTyp    
  end interface absTyp

contains

!!--------------------------initAbsTyp-------------------------!!
  function initAbsTyp(N,typ,x0,l,wT) result(b)
  implicit none
    
    integer(kind=C_K1),intent(in)::N,typ
    real(kind=C_K2),intent(in)::x0,l,wT    
    type(absTyp)::b    

    b%N=N
    b%typ=typ
    b%x0=x0
    b%l=l
    b%wT=wT
    b%c1=30d0/wT/(dexp(1d0)-1d0)    

  end function initAbsTyp
!!------------------------End initAbsTyp-----------------------!!



!!---------------------------calcAbsC--------------------------!!
  subroutine calcAbsC(b,npt,cor,absC)
  implicit none

    class(absTyp),intent(in)::b
    integer(kind=C_K1),intent(in)::npt
    integer(kind=C_K1)::i,ix

    real(kind=C_K2),intent(in)::cor(npt,2)
    real(kind=C_K2),intent(inout)::absC(npt)
    real(kind=C_K2)::x,dc,dx

    select case(b%typ)
      case(1)
        ix=2
        dc=1d0

      case(2)
        ix=1
        dc=-1d0

      case(3)
        ix=2
        dc=-1d0

      case(4)
        ix=1
        dc=1d0

    end select

    do i=1,npt
      dx=dc*(cor(i,ix)-b%x0)/b%l
      if(dx.gt.0d0)then
        absC(i)=b%c1*(dexp(dx**2)-1d0)
      endif
    enddo

  end subroutine calcAbsC
!!-------------------------End calcAbsC------------------------!!

end module outAbsModule
!!------------------------End outAbsModule-------------------------!!



!!-------------------------waveFileModule--------------------------!!
module waveFileModule
use bsnqGlobVars
implicit none

  type, public :: wvFileType        
    character(len=256)::fileName
    integer(kind=C_K1)::numP,posI
    real(kind=C_K2),allocatable::data(:,:)
  contains
    procedure ::  initWaveFile
    procedure ::  getEta
    procedure ::  getPQ
  end type wvFileType  

contains

!!------------------------initWaveFile-------------------------!!
  subroutine initWaveFile(b)
  implicit none

    class(wvFileType),intent(inout)::b
    integer(kind=C_K1)::mf,i
    real(kind=C_K2)::tmpra(4)
    logical::ex

    inquire(file=trim(b%fileName),exist=ex)
    if(ex) then
      open(newunit=mf,file=trim(b%fileName))
    else
      write(9,*)"[ERR] Missing wave input file"
      stop
    endif

    b%numP=0
    do while(.true.)
      read(mf,*,end=11,err=12)tmpra
      b%numP=b%numP+1
      cycle
      11 exit
      12 write(9,*)"[ERR] Check wave file format"
      stop
    enddo
    write(9,'(" [INF] Number of wave file point = ",i10)')b%numP
    close(mf)

    allocate(b%data(b%numP,4))
    open(newunit=mf,file=trim(b%fileName))    
    do i=1,b%numP
      read(mf,*,end=21,err=21)b%data(i,1:4)      
    enddo
    write(9,*)"[INF] Done wave file read"

    b%posI=2

    return
    21 write(9,*)"[ERR] Check wave file format"
    stop
  end subroutine initWaveFile
!!----------------------End initWaveFile-----------------------!!




!!---------------------------getEta----------------------------!!
  subroutine getEta(b,rTime,eta)
  implicit none

    class(wvFileType),intent(inout)::b

    integer(kind=C_K1)::k

    real(kind=C_K2),intent(in)::rTime
    real(kind=C_K2),intent(out)::eta
    
    if((rTime.gt.b%data(b%numP,1)).or.&
      (rTime.lt.b%data(1,1)))then
      write(9,*)"[ERR] Wave query exceeds supplied time range"
      write(9,'( "[---] ",F15.6)')rTime
      write(9,'( "[---] ",2F15.6)')b%data(1,1),b%data(b%numP,1)
      stop
    endif

    do while(b%data(b%posI+1,1).lt.rTime)
      b%posI=b%posI+1
    enddo

    do while(b%data(b%posI-1,1).gt.rTime)
      if(b%posI.eq.2)exit
      b%posI=b%posI-1
    enddo

    k=2
    eta=b%data(b%posI,k) &
      + (b%data(b%posI+1,k)-b%data(b%posI-1,k)) &
      / (b%data(b%posI+1,1)-b%data(b%posI-1,1)) &
      * (rTime-b%data(b%posI,1))

  end subroutine getEta
!!-------------------------End getEta--------------------------!!



!!----------------------------getPQ----------------------------!!
  subroutine getPQ(b,rTime,p,q)
  implicit none

    class(wvFileType),intent(inout)::b

    integer(kind=C_K1)::k

    real(kind=C_K2),intent(in)::rTime
    real(kind=C_K2),intent(out)::p,q
    
    if((rTime.gt.b%data(b%numP,1)).or.&
      (rTime.lt.b%data(1,1)))then
      write(9,*)"[ERR] Wave query exceeds supplied time range"
      write(9,'( "[---] ",F15.6)')rTime
      write(9,'( "[---] ",2F15.6)')b%data(1,1),b%data(b%numP,1)
      stop
    endif

    do while(b%data(b%posI+1,1).lt.rTime)
      b%posI=b%posI+1
    enddo

    do while(b%data(b%posI-1,1).gt.rTime)
      if(b%posI.eq.2)exit
      b%posI=b%posI-1
    enddo

    k=3
    p=b%data(b%posI,k) &
      + (b%data(b%posI+1,k)-b%data(b%posI-1,k)) &
      / (b%data(b%posI+1,1)-b%data(b%posI-1,1)) &
      * (rTime-b%data(b%posI,1))

    k=4
    q=b%data(b%posI,k) &
      + (b%data(b%posI+1,k)-b%data(b%posI-1,k)) &
      / (b%data(b%posI+1,1)-b%data(b%posI-1,1)) &
      * (rTime-b%data(b%posI,1))

  end subroutine getPQ
!!--------------------------End getPQ--------------------------!!


end module waveFileModule
!!-----------------------End waveFileModule------------------------!!