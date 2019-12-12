!!-----------------------------shipMod-----------------------------!!
module shipMod
use bsnqGlobVars
implicit none

  type, public :: shipType    
    integer(kind=C_K1)::posN,posI,totNShip
    real(kind=C_K2)::cl,cb,al    
    real(kind=C_K2)::L,B,T
    real(kind=C_K2),allocatable::posData(:,:)
  contains
    procedure ::  getPress    
    procedure ::  getLoc
  end type shipType

  interface shipType
    procedure :: initShip
  end interface shipType
  
contains



!!------------------------initShipType-------------------------!!
  type(shipType) function initShip(mf,numShip,simEnd)
  implicit none

    integer(kind=C_K1),intent(in)::mf,numShip
    real(kind=C_K2),intent(in)::simEnd
    integer(kind=C_K1)::i,j
    character(len=C_KSTR)::bqtxt  

    initShip%totNShip=numShip
    read(mf,*,end=83,err=83)bqtxt
    read(mf,*,end=83,err=83)bqtxt
    read(mf,*,end=83,err=83)initShip%cl,initShip%cb,initShip%al
    read(mf,*,end=83,err=83)bqtxt
    read(mf,*,end=83,err=83)initShip%L,initShip%B,initShip%T
    read(mf,*,end=83,err=83)bqtxt
    read(mf,*,end=83,err=83)initShip%posN        
    read(mf,*,end=83,err=83)bqtxt
    allocate(initShip%posData( initShip%posN, 4))
    do j=1,initShip%posN
      read(mf,*,end=83,err=84)initShip%posData(j,1:4)
    enddo        
    initShip%posI=1

    if(simEnd.gt.initShip%posData( initShip%posN,1 ))then
      write(9,*)"[ERR] Insufficient ship position time information"
      stop
    endif    

    goto 84
    83 write(9,*) "[ERR] Check ship file format"
    stop
    84 continue

  end function initShip
!!----------------------End initShipType-----------------------!!




!!---------------------------getLoc----------------------------!!
  subroutine getLoc(sh,rTime,x0,y0,thDeg)
  implicit none

    class(shipType),intent(inout)::sh    
    integer(kind=C_K1)::i,k,lPosI
    real(kind=C_K2),intent(in)::rTime
    real(kind=C_K2),intent(out)::x0,y0,thDeg

    if((rTime.ge.sh%posData(sh%posI,1)) .and. &
      (rTime.lt.sh%posData(sh%posI+1,1)))then
      lPosI=sh%posI      
    
    else
      do i=1,sh%posN
        if(sh%posData(i,1).gt.rTime) exit      
      enddo
      lPosI=i-1
      if((lPosI.eq.0) .or. (rTime.gt.sh%posData(sh%posN,1)) )then
        write(9,*)"[ERR] Ship position unavailable at time ",rTime
        stop
      endif
      sh%posI=lPosI
    endif

    k=2
    x0=sh%posData(lPosI,k) &
      +(sh%posData(lPosI+1,k)-sh%posData(lPosI,k)) &
      /(sh%posData(lPosI+1,1)-sh%posData(lPosI,1)) &
      *(rTime-sh%posData(lPosI,1))
    k=3
    y0=sh%posData(lPosI,k) &
      +(sh%posData(lPosI+1,k)-sh%posData(lPosI,k)) &
      /(sh%posData(lPosI+1,1)-sh%posData(lPosI,1)) &
      *(rTime-sh%posData(lPosI,1))
    k=4
    thDeg=sh%posData(lPosI,k) &
      +(sh%posData(lPosI+1,k)-sh%posData(lPosI,k)) &
      /(sh%posData(lPosI+1,1)-sh%posData(lPosI,1)) &
      *(rTime-sh%posData(lPosI,1))

  end subroutine getLoc
!!-------------------------End getLoc--------------------------!!




!!--------------------------getPress---------------------------!!
  subroutine getPress(sh,rTime,npt,cor,pres)
  implicit none

    class(shipType),intent(inout)::sh
    integer(kind=C_K1),intent(in)::npt
    real(kind=C_K2),intent(in)::rTime,cor(npt,2)
    real(kind=C_K2),intent(out)::pres(npt)
    
    integer(kind=C_K1)::i
    real(kind=C_K2)::x0,y0,thDeg,thRad,rRef,cost,sint
    real(kind=C_K2)::x,y,dr2,lx,ly
    
    call sh%getLoc(rTime,x0,y0,thDeg)
    thRad=thDeg*deg2rad

    pres=0d0
    rRef=(max(sh%L,sh%B)/2d0*1.2d0)**2
    cost=dcos(thRad)
    sint=dsin(thRad)

    do i=1,npt
      x=cor(i,1)-x0
      y=cor(i,2)-y0
      dr2=(x**2 + y**2)
      if(dr2.lt.rRef)then
        lx=(+x*cost + y*sint)/sh%L
        ly=(-x*sint + y*cost)/sh%B
        ! lx=(-x*cost - y*sint)/L
        ! ly=(+x*sint - y*cost)/B
        if(abs(lx).gt.0.5d0)cycle
        if(abs(ly).gt.0.5d0)cycle        
        pres(i)=sh%T*(1-sh%cl*(lx**4))*(1-sh%cb*(ly**2)) &
          *exp(-sh%al*(ly**2))
      endif
    enddo

  end subroutine getPress
!!------------------------End getPress-------------------------!!


end module shipMod
!!---------------------------End shipMod---------------------------!!