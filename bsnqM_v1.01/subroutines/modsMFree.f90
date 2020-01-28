!!---------------------------meshFreeMod---------------------------!!
module meshFreeMod
use bsnqGlobVars
implicit none
  
  type, public :: mfTyp
    real(kind=C_K2),allocatable::rad(:),phiDx(:),phiDy(:)
  contains
    !procedure :: initRadius
    procedure :: calcRadius
  end type mfTyp

contains


! !!-------------------------initRadius--------------------------!!
!   subroutine initRadius(m,np)
!   implicit none

!     class(mfTyp),intent(inout)::m
!     integer(kind=C_K1),intent(in)::np

!     if(allocated(m%rad)) deallocate(m%rad)
!     allocate(m%rad(np))
!     m%rad=0d0

!   end subroutine initRadius
! !!-----------------------End initRadius------------------------!!



!!-------------------------calcRadius--------------------------!!
  subroutine calcRadius(m,npt,sz,ivq,jvq,corx,cory)
  implicit none

    class(mfTyp),intent(inout)::m
    integer(kind=C_K1),intent(in)::npt,sz
    integer(kind=C_K1),intent(in)::ivq(0:npt),jvq(sz)
    real(kind=C_K2),intent(in)::corx(npt),cory(npt)

    integer(kind=C_K1)::i,j,k1,k2,i2
    real(kind=C_K2)::rMax,dr,cx,cy,coef

    coef=0.8d0 ! Farhter point in this much of the radius

    if(allocated(m%rad)) deallocate(m%rad)
    allocate(m%rad(npt))    
    m%rad=0d0

    if(allocated(m%phiDx)) deallocate(m%phiDx)
    if(allocated(m%phiDy)) deallocate(m%phiDy)
    allocate(m%phiDx(sz),m%phiDy(sz))
    m%phiDx=0d0
    m%phiDy=0d0
    
    ! Radius
    do i=1,npt
      rMax=0d0
      cx=corx(i)
      cy=cory(i)
      k1=(i-1)*ivq(0)
      k2=k1+ivq(i)
      do j=k1+1,k2
        i2=jvq(j)
        dr=(corx(i2)-cx)**2 + (cory(i2)-cy)**2
        if(rMax.lt.dr) rMax=dr
      enddo
      if(rMax.lt.1e-10)then
        write(9,'(" [ERR] Check radius calculation at node",I10)')i
        stop
      endif
      m%rad(i)=dsqrt(rMax)/coef
    enddo

    do i=1,npt
      k1=(i-1)*ivq(0)+1
      k2=(i-1)*ivq(0)+ivq(i)-1
      call sfdi2DDx(corx(i),cory(i),ivq(i)-1,m%rad(i),&
        corx(jvq(k1:k2)),cory(jvq(k1:k2)),&
        m%phiDx(k1:k2),m%phiDy(k1:k2))      
      stop
    enddo

  end subroutine calcRadius
!!-----------------------End calcRadius------------------------!!



!!--------------------------sfdi2DDx---------------------------!!
  subroutine sfdi2DDx(xi,yi,nn,riav,corx,cory,phiDx,phiDy)
  implicit none

    integer(kind=C_K1),intent(in)::nn    
    real(kind=C_K2),intent(in)::xi,yi,corx(nn),cory(nn),riav
    real(kind=C_K2),intent(out)::phiDx(nn),phiDy(nn)

    integer(kind=C_K1)::i,j,k1,k2,i2
    real(kind=C_K2)::rMax,dr,dx,dy,wj,drr
    real(kind=C_K2)::rnxx,rnxy,rnyy
    real(kind=C_K2)::bjx,bjy

    phiDx=0d0
    phiDy=0d0

    rnxx=0d0
    rnyy=0d0
    rnxy=0d0

    do i=1,nn
      dx=corx(i)-xi
      dy=cory(i)-yi
      dr=dsqrt(dx**2 + dy**2)
      drr=dr/riav

      wj=0d0
      if((drr.gt.1e-10).and.(drr.le.1d0))then
        wj=1d0 - 6d0*drr**2 + 8*drr**3 - 3*drr**4
      endif

      if(wj.gt.1e-15)then
        rnxx=rnxx + dx*dx*wj/dr/dr
        rnxy=rnxy + dx*dy*wj/dr/dr
        rnyy=rnyy + dy*dy*wj/dr/dr
      endif
    enddo

    do i=1,nn
      dx=corx(i)-xi
      dy=cory(i)-yi
      dr=dsqrt(dx**2 + dy**2)
      drr=dr/riav

      wj=0d0
      if((drr.gt.1e-10).and.(drr.le.1d0))then
        wj=1d0 - 6d0*drr**2 + 8*drr**3 - 3*drr**4
      endif

      if(wj.gt.1e-15)then
        bjx=-wj*dx/dr/dr
        bjy=-wj*dy/dr/dr

        phiDx(i)=bjx/rnxx
        phiDx(i)=bjx/rnxx        
      endif
    enddo

  end subroutine sfdi2DDx
!!------------------------End sfdi2DDx-------------------------!!


end module meshFreeMod
!!-------------------------End meshFreeMod-------------------------!!



