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
  
contains

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
    
    real(kind=C_K2)::x0,y0,thDeg,thRad
    
    call sh%getLoc(rTime,x0,y0,thDeg)
    thRad=thDeg*deg2rad




  end subroutine getPress
!!------------------------End getPress-------------------------!!


end module shipMod
!!---------------------------End shipMod---------------------------!!