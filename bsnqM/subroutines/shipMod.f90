!!-----------------------------shipMod-----------------------------!!
module shipMod
use bsnqGlobVars
use meshFreeMod
implicit none

  type, public :: shipType    
    integer(kind=C_K1)::posN,posI,totNShip,presFnc
    integer(kind=C_K1)::gridNL,gridNB
    real(kind=C_K2)::cl,cb,al    
    real(kind=C_K2)::L,B,T
    real(kind=C_K2),allocatable::posData(:,:)
    real(kind=C_K2),allocatable::gP(:,:) !X major form
    real(kind=C_K2),allocatable::gPPres(:)
    logical(kind=C_LG)::dragFlag, initEtaFlag
    type(mfPoiTyp)::gPObj    
  contains
    procedure ::  getPress    
    procedure ::  getLoc
    procedure ::  calcDrag
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
    integer(kind=C_K1)::i,j,nnMax,k,k2
    real(kind=C_K2)::tmpr1,tmpr2
    character(len=C_KSTR)::bqtxt  

    initShip%totNShip=numShip
    read(mf,*,end=83,err=83)bqtxt
    read(mf,*,end=83,err=83)bqtxt
    read(mf,*,end=83,err=83)initShip%presFnc
    select case(initShip%presFnc)
      case (1)        
        read(mf,*,end=83,err=83)bqtxt
        read(mf,*,end=83,err=83)initShip%cl,initShip%cb

      case (2) 
        read(mf,*,end=83,err=83)bqtxt
        read(mf,*,end=83,err=83)initShip%cl,initShip%cb,initShip%al
      
      case (3)
        initShip%cl=0
        initShip%cb=0
        initShip%al=0
      
      case DEFAULT
        write(9,*)'[ERR] Invalid ship type'
        stop
    end select

    read(mf,*,end=83,err=83)bqtxt
    read(mf,*,end=83,err=83)initShip%L,initShip%B,initShip%T

    read(mf,*,end=83,err=83)bqtxt
    read(mf,*,end=83,err=83)initShip%initEtaFlag

    read(mf,*,end=83,err=83)bqtxt
    read(mf,*,end=83,err=83)bqtxt
    read(mf,*,end=83,err=83)initShip%dragFlag
    read(mf,*,end=83,err=83)bqtxt
    read(mf,*,end=83,err=83)nnMax,tmpr1
    read(mf,*,end=83,err=83)bqtxt
    read(mf,*,end=83,err=83)initShip%gridNL,initShip%gridNB
    allocate( initShip%gP( initShip%gridNL*initShip%gridNB, 2 ) )
    allocate( initShip%gPPres( initShip%gridNL*initShip%gridNB ) )
    call initShip%gPObj%initPoi(nnMax,0d0,0d0,tmpr1)    

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
    real(kind=C_K2)::x,y,dr2,lx,ly,px,py
    
    call sh%getLoc(rTime,x0,y0,thDeg)
    thRad=thDeg*deg2rad

    pres=0d0

    select case(sh%presFnc)
    case (1)
        rRef=(max(sh%L,sh%B)/2d0*1.5d0)**2
        cost=dcos(thRad)
        sint=dsin(thRad)

        do i=1,npt
          x=cor(i,1)-x0
          y=cor(i,2)-y0
          dr2=(x**2 + y**2)
          if(dr2.lt.rRef)then
            lx=abs(+x*cost + y*sint)/sh%L
            ly=abs(-x*sint + y*cost)/sh%B
            ! lx=(-x*cost - y*sint)/L
            ! ly=(+x*sint - y*cost)/B
            if(lx.gt.0.5d0)cycle
            if(ly.gt.0.5d0)cycle
            
            if(lx.lt.0.5d0*sh%cl)then
              px=1        
            else
              px=dcos(pi*(lx-0.5d0*sh%cl)/(1d0-sh%cl))      
            endif
            
            if(ly.lt.0.5d0*sh%cb)then
              py=1
            else
              py=dcos(pi*(ly-0.5d0*sh%cb)/(1d0-sh%cb))      
            endif
            
            pres(i)=sh%T*px*px*py*py
          endif
        enddo

      case (2)
        rRef=(max(sh%L,sh%B)/2d0*1.5d0)**2
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

      case (3)        
        do i=1,npt
          x=abs(cor(i,1)-x0)/sh%L
          if(x.gt.0.5d0)cycle
          lx=dcos(pi*x)
          pres(i)=sh%T*lx**2
        enddo

    end select

  end subroutine getPress
!!------------------------End getPress-------------------------!!



!!--------------------------calcDrag---------------------------!!
  subroutine calcDrag(sh,rTime,np,corx,cory,eta,fx)
  implicit none

    class(shipType),intent(inout)::sh
    integer(kind=C_K1),intent(in)::np
    real(kind=C_K2),intent(in)::rTime,corx(np),cory(np)
    real(kind=C_K2),intent(in)::eta(np)
    real(kind=C_K2),intent(out)::fx

    integer(kind=C_K1)::i,j,k,k2,shI,shJ,nn,err
    real(kind=C_K2)::x0,y0,thDeg,csth,snth,dr,tmpr1,etaDx
    real(kind=C_K2)::dx,dy,x,y

    call sh%getLoc(rTime,x0,y0,thDeg)
    csth=dcos(thDeg*deg2rad)
    snth=dsin(thDeg*deg2rad)

    dx = sh%L / (sh%gridNL-1)
    dy = sh%B / (sh%gridNB-1)
    k2=sh%gridNB
    do i=1,sh%gridNL
      do j=1,sh%gridNB
        k=(i-1)*k2+j
        x= -sh%L/2d0 + (i-1)*dx
        y= -sh%B/2d0 + (j-1)*dy
        sh%gP(k,1) = x0 + x*csth - y*snth  
        sh%gP(k,2) = y0 + x*snth + y*csth  
      enddo
    enddo
    ! do j=sh%gridNB,1,-1
    !   do i=1,sh%gridNL
    !     k=(i-1)*k2+j
    !     write(*,'(F6.2,",",F6.2," ")',advance='no')sh%gP(k,1), &
    !       sh%gP(k,2)
    !   enddo
    !   write(*,*)
    ! enddo
    ! stop

    k2=sh%gridNL*sh%gridNB
    call sh%getPress(rTime,k2,sh%gP,sh%gPPres)

    fx=0d0
    k2=sh%gridNB
    do shI=1,sh%gridNL
      do shJ=1,sh%gridNB

        k=(shI-1)*k2+shJ
        sh%gPObj%cx = sh%gP(k,1)
        sh%gPObj%cy = sh%gP(k,2)

        nn=0
        tmpr1=sh%gPObj%rad**2
        do i=1,np
          dr=( corx(i)-sh%gPObj%cx )**2 &
            + ( cory(i)-sh%gPObj%cy )**2 
          if(dr.lt.tmpr1)then
            nn=nn+1
            if(nn .gt. sh%gPObj%nnMax)then
              write(9,'("      |",a6,a)')"[ERR]",&
                "ship%CalcDrag| Increase ship numOfNegh or reduce radius"
              fx=0d0
              return
            endif            
            sh%gPObj%neid(nn) = i
          endif
        enddo
        if(nn .lt. 4)then
          write(9,'("      |",a6,a)')"[ERR]",&
            "ship%CalcDrag| Insufficient neghs"
          ! write(*,'(2I5,2F15.6)')shI,shJ,sh%gPObj%cx,sh%gPObj%cy
          ! write(*,'(I5)')nn
          fx=0d0
          return
        endif            
        sh%gPObj%nn=nn        

        call mls2DDx(sh%gPObj%cx, sh%gPObj%cy, nn, &
          sh%gPObj%rad, corx(sh%gPObj%neid(1:nn)), &
          cory(sh%gPObj%neid(1:nn)), &
          sh%gPObj%phi(1:nn), sh%gPObj%phiDx(1:nn), &
          sh%gPObj%phiDy(1:nn),err)        

        if(err.ne.0) then !If any err in mls2DDx then return fx=0
          fx=0d0
          return
        endif

        !! etaDx
        etaDx=0d0
        do i=1,nn
          etaDx=etaDx+ sh%gPObj%phiDx(i) * eta(sh%gPObj%neid(i))
        enddo

        !! Simpson's integration
        tmpr1=1d0
        if(mod(shI,2).ne.0)then
          if((shI.eq.1).or.(shI.eq.sh%gridNL))then
            tmpr1=tmpr1*1d0
          else
            tmpr1=tmpr1*2d0
          endif
        else
          tmpr1=tmpr1*4d0
        endif
        if(mod(shJ,2).ne.0)then
          if((shJ.eq.1).or.(shJ.eq.sh%gridNB))then
            tmpr1=tmpr1*1d0
          else
            tmpr1=tmpr1*2d0
          endif
        else
          tmpr1=tmpr1*4d0
        endif
        tmpr1=tmpr1*dx*dy/9d0
        fx=fx+tmpr1*sh%gPPres(k)*etaDx


      enddo
    enddo

  end subroutine calcDrag
!!------------------------End calcDrag-------------------------!!


end module shipMod
!!---------------------------End shipMod---------------------------!!