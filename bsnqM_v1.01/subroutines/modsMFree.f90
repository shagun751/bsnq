!!---------------------------meshFreeMod---------------------------!!
module meshFreeMod
use bsnqGlobVars
implicit none
  
  type, public :: mfTyp
    real(kind=C_K2),allocatable::rad(:),phi(:)
    real(kind=C_K2),allocatable::phiDx(:),phiDy(:)
  contains
    !procedure :: initRadius
    procedure :: calcAll
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



!!---------------------------calcAll---------------------------!!
  subroutine calcAll(m,npt,sz,ivq,jvq,corx,cory)
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

    if(allocated(m%phi)) deallocate(m%phi)
    if(allocated(m%phiDx)) deallocate(m%phiDx)
    if(allocated(m%phiDy)) deallocate(m%phiDy)
    allocate(m%phi(sz), m%phiDx(sz), m%phiDy(sz))
    m%phi=0d0
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
      k2=(i-1)*ivq(0)+ivq(i)
      call mls2DDx(corx(i), cory(i), ivq(i), m%rad(i), &
        corx(jvq(k1:k2)), cory(jvq(k1:k2)), &
        m%phi(k1:k2), m%phiDx(k1:k2), m%phiDy(k1:k2))      

      ! write(*,*)ivq(i)
      ! write(*,*)corx(i),cory(i)
      ! write(*,*)corx(jvq(k1:k2))
      ! write(*,*)cory(jvq(k1:k2))
      ! write(*,*)(m%phi(k1:k2))
      ! write(*,*)(m%phiDx(k1:k2))
      ! write(*,*)(m%phiDy(k1:k2))
      ! stop
    enddo
    

  end subroutine calcAll
!!-------------------------End calcAll-------------------------!!



!!-------------------------testMls2DDx-------------------------!!
  subroutine testMls2DDx
  implicit none

    integer(kind=C_K1)::nn,i
    real(kind=C_K2)::xi,yi,rad,tmpr1,tmpr2,tmpr3
    real(kind=C_K2),allocatable::corx(:),cory(:)
    real(kind=C_K2),allocatable::phi(:),phiDx(:),phiDy(:)
    real(kind=C_K2)::A(3,3),AInv(3,3),ADet

    xi=0d0
    yi=0d0
    nn=9
    rad=dsqrt(2d0)/0.8d0

    allocate(corx(nn),cory(nn),phi(nn),phiDx(nn),phiDy(nn))

    corx = (/ 1d0, 1d0, 0d0, -1d0, -1d0, -1d0, 0d0, 1d0, 0d0 /)
    cory = (/ 0d0, 1d0, 1d0, 1d0, 0d0, -1d0, -1d0, -1d0, 0d0 /)

    call mls2DDx(xi, yi, nn, rad, &
        corx, cory, &
        phi, phiDx, phiDy)          

    write(*,'(4F15.6)')phi
    write(*,'(4F15.6)')phiDx
    write(*,'(4F15.6)')phiDy
    write(*,'(3F15.6)')sum(phi),sum(phiDx),sum(phiDy)

    tmpr1=0d0
    tmpr2=0d0
    tmpr3=0d0
    do i=1,nn
      tmpr1 = tmpr1 + phi(i)*(1+corx(i))*(10+cory(i))
    enddo
    write(*,'(3F15.6)')tmpr1,tmpr2,tmpr3

    ! A(1,:) = (/ 1d0,0.2d0,3d0 /)
    ! A(2,:) = (/ 0.2d0,4d0,5d0 /)
    ! A(3,:) = (/ 3d0,5d0,6d0 /)

    ! call findInvSymm3x3(A,AInv,ADet)
    ! write(*,*)
    ! write(*,'(3F15.6)')AInv
    ! write(*,'(F15.6)')ADet

  end subroutine testMls2DDx
!!-----------------------End testMls2DDx-----------------------!!



!!---------------------------mls2DDx---------------------------!!
  subroutine mls2DDx(x,y,nn,R,corx,cory,phi,phiDx,phiDy)
  implicit none

    integer(kind=C_K1),intent(in)::nn    
    real(kind=C_K2),intent(in)::x,y,corx(nn),cory(nn),R
    real(kind=C_K2),intent(out)::phi(nn),phiDx(nn),phiDy(nn)

    integer(kind=C_K1)::j,k1,k2,i2
    real(kind=C_K2)::rMax,dr,dx,dy,wj,wjx,wjy,drr,xj,yj
    real(kind=C_K2)::A(3,3),AInv(3,3),ADet
    real(kind=C_K2)::V13(3),V31(3),M33(3,3)

    phi=0d0
    phiDx=0d0
    phiDy=0d0
   
    A=0d0    
    do j=1,nn
      xj=corx(j)
      yj=cory(j)
      dx=xj-x
      dy=yj-y
      dr=dsqrt(dx**2 + dy**2)
      drr=dr/R

      if(drr.le.1d0)then
        ! !! Exp (ending at r/R = 1)
        ! wj=0.5d0/pi*dexp(-4.5d0*drr*drr)
        ! wjx=9d0*dx/R/R*wj
        ! wjy=9d0*dy/R/R*wj

        !! Biquadratic
        wj = 1d0 - 6d0*drr**2 + 8*drr**3 - 3*drr**4        
      endif

      A(1,1)=A(1,1) + wj
      A(1,2)=A(1,2) + wj*xj
      A(1,3)=A(1,3) + wj*yj
      A(2,2)=A(2,2) + wj*xj*xj
      A(2,3)=A(2,3) + wj*xj*yj
      A(3,3)=A(3,3) + wj*yj*yj        
    enddo    
    A(2,1)=A(1,2)
    A(3,1)=A(1,3)
    A(3,2)=A(2,3)
    call findInvSymm3x3(A,AInv,ADet)

    if(abs(ADet).lt.1e-10)then
      write(9,'(" [ERR] ADet is too small",F15.6)')ADet
      write(9,'(" [---] Location X Y",2F15.6)')x,y
      write(9,'(" [---] Num Neigh",I10)')nn
      stop
    endif    

    do j=1,nn
      xj=corx(j)
      yj=cory(j)
      dx=xj-x
      dy=yj-y
      dr=dsqrt(dx**2 + dy**2)
      drr=dr/R

      if(drr.le.1d0)then
        ! !! Exp (ending at r/R = 1)
        ! wj=0.5d0/pi*dexp(-4.5d0*drr*drr)
        ! wjx=9d0*dx/R/R*wj
        ! wjy=9d0*dy/R/R*wj

        !! Biquadratic
        wj = 1d0 - 6d0*drr**2 + 8*drr**3 - 3*drr**4        
      endif


      V13(1)=1d0
      V13(2)=x
      V13(3)=y

      V31(1)=wj
      V31(2)=wj*xj
      V31(3)=wj*yj

      call matMul_V13_M33_V31(V13,AInv,V31,phi(j))
            
    enddo    

  end subroutine mls2DDx
!!------------------------End mls2DDx-------------------------!!



!!-----------------------mls2DDxSAThesis-----------------------!!
  subroutine mls2DDxSAThesis(xi,yi,nn,R,corx,cory,phi,phiDx,phiDy)
  implicit none

    integer(kind=C_K1),intent(in)::nn    
    real(kind=C_K2),intent(in)::xi,yi,corx(nn),cory(nn),R
    real(kind=C_K2),intent(out)::phi(nn),phiDx(nn),phiDy(nn)

    integer(kind=C_K1)::j,k1,k2,i2
    real(kind=C_K2)::rMax,dr,dx,dy,wj,wjx,wjy,drr
    real(kind=C_K2)::m(0:5),mx(0:5),my(0:5)
    real(kind=C_K2)::c0,c1,c21,c22
    real(kind=C_K2)::c0x,c1x,c21x,c22x
    real(kind=C_K2)::c0y,c1y,c21y,c22y

    !m(0:5) = (/ m0 m1 m2 m11 m12 m22 /)

    phi=0d0
    phiDx=0d0
    phiDy=0d0
   
    m=0d0
    mx=0d0
    my=0d0
    do j=1,nn
      dx=corx(j)-xi
      dy=cory(j)-yi
      dr=dsqrt(dx**2 + dy**2)
      drr=dr/R

      wj=0d0      
      if(drr.le.1d0)then
        !! Exp (ending at r/R = 1)
        wj=0.5d0/pi*dexp(-4.5d0*drr*drr)
        wjx=9d0*dx/R/R*wj
        wjy=9d0*dy/R/R*wj

        !! Biquadratic
        ! wj = 1d0 - 6d0*drr**2 + 8*drr**3 - 3*drr**4
        ! wjx= 12d0*dx/(R**2) - 24d0*dx*dr/(R**3) &
        !   + 12d0*dx*dr*dr/(R**4)
        ! wjy= 12d0*dy/(R**2) - 24d0*dy*dr/(R**3) &
        !   + 12d0*dy*dr*dr/(R**4)
      endif

      !if(wj.gt.1e-15)then
        m(0)=m(0)+wj
        m(1)=m(1)+wj*dx/R
        m(2)=m(2)+wj*dy/R
        m(3)=m(3)+wj*dx*dx/R/R
        m(4)=m(4)+wj*dx*dy/R/R
        m(5)=m(5)+wj*dy*dy/R/R

        mx(0)=mx(0)+wjx
        mx(1)=mx(1)+wjx*dx/R - wj/R
        mx(2)=mx(2)+wjx*dy/R
        mx(3)=mx(3)+wjx*dx*dx/R/R - 2d0*wj*dx/R/R
        mx(4)=mx(4)+wjx*dx*dy/R/R - wj*dy/R/R
        mx(5)=mx(5)+wjx*dy*dy/R/R

        my(0)=my(0)+wjy
        my(1)=my(1)+wjy*dx/R
        my(2)=my(2)+wjy*dy/R - wj/R
        my(3)=my(3)+wjy*dx*dx/R/R
        my(4)=my(4)+wjy*dx*dy/R/R - wj*dx/R/R
        my(5)=my(5)+wjy*dy*dy/R/R - 2d0*wj*dy/R/R
      !endif
    enddo    

    c0 = m(0)*( m(4)*m(6) - m(5)**2 ) - ( m(1)*m(1)*m(6) &
      - 2d0*m(1)*m(2)*m(5) + m(2)*m(2)*m(4) )

    if(abs(c0).lt.1e-10)then
      write(9,'(" [ERR] C0 is too small",F15.6)')c0
      write(9,'(" [---] Location X Y",2F15.6)')xi,yi
      write(9,'(" [---] Num Neigh",I10)')nn
      stop
    endif

    c1 = (m(4)*m(6) - m(5)**2)/c0
    c21 = (m(2)*m(5) - m(1)*m(6))/c0
    c22 = (m(1)*m(5) - m(2)*m(4))/c0

    c0x = m(0)*( mx(4)*m(6) + m(4)*mx(6) - 2d0*m(5)*mx(5) ) &
      + mx(0)*( m(4)*m(6) - m(5)**2 ) &
      - ( 2d0*m(1)*mx(1)*m(6) + m(1)*m(1)*mx(6) ) &
      - ( 2d0*m(2)*mx(2)*m(4) + m(2)*m(2)*mx(4) ) &
      + ( 2d0*mx(1)*m(2)*m(5) + 2d0*m(1)*mx(2)*m(5) &
        + 2d0*m(1)*m(2)*mx(5) )
    c1x = ( mx(4)*m(6) + m(4)*mx(6) - 2d0*mx(5)*m(5) )/c0 &
      - ( m(4)*m(6) - m(5)**2 )*c0x/c0/c0
    c21x = ( mx(2)*m(5) + m(2)*mx(5) - mx(1)*m(6) - m(1)*mx(6) )/c0 &
      - ( m(2)*m(5) - m(1)*m(6) )*c0x/c0/c0
    c22x = ( mx(1)*m(5) + m(1)*mx(5) - mx(2)*m(4) - m(2)*mx(4) )/c0 &
      - ( m(1)*m(5) - m(2)*m(4) )*c0x/c0/c0

    c0y = m(0)*( my(4)*m(6) + m(4)*my(6) - 2d0*m(5)*my(5) ) &
      + my(0)*( m(4)*m(6) - m(5)**2 ) &
      - ( 2d0*m(1)*my(1)*m(6) + m(1)*m(1)*my(6) ) &
      - ( 2d0*m(2)*my(2)*m(4) + m(2)*m(2)*my(4) ) &
      + ( 2d0*my(1)*m(2)*m(5) + 2d0*m(1)*my(2)*m(5) &
        + 2d0*m(1)*m(2)*my(5) )
    c1y = ( my(4)*m(6) + m(4)*my(6) - 2d0*my(5)*m(5) )/c0 &
      - ( m(4)*m(6) - m(5)**2 )*c0y/c0/c0
    c21y = ( my(2)*m(5) + m(2)*my(5) - my(1)*m(6) - m(1)*my(6) )/c0 &
      - ( m(2)*m(5) - m(1)*m(6) )*c0y/c0/c0
    c22y = ( my(1)*m(5) + m(1)*my(5) - my(2)*m(4) - m(2)*my(4) )/c0 &
      - ( m(1)*m(5) - m(2)*m(4) )*c0y/c0/c0


    do j=1,nn
      dx=corx(j)-xi
      dy=cory(j)-yi
      dr=dsqrt(dx**2 + dy**2)
      drr=dr/R

      wj=0d0      
      if(drr.le.1d0)then
        !! Exp (ending at r/R = 1)
        wj=0.5d0/pi*dexp(-4.5d0*drr*drr)
        wjx=9d0*dx/R/R*wj
        wjy=9d0*dy/R/R*wj

        !! Biquadratic
        ! wj = 1d0 - 6d0*drr**2 + 8*drr**3 - 3*drr**4
        ! wjx= 12d0*dx/(R**2) - 24d0*dx*dr/(R**3) &
        !   + 12d0*dx*dr*dr/(R**4)
        ! wjy= 12d0*dy/(R**2) - 24d0*dy*dr/(R**3) &
        !   + 12d0*dy*dr*dr/(R**4)
      endif

      !if(wj.gt.1e-15)then
        phi(j)=wj*( c1 + c21*dx/R + c22*dy/R )
        phiDx(j)=wj*( c1x + c21x*dx/R + c22x*dy/R -c21/R ) &
          + wjx*( c1 + c21*dx/R + c22*dy/R )
        phiDy(j)=wj*( c1y + c21y*dx/R + c22y*dy/R -c22/R ) &
          + wjy*( c1 + c21*dx/R + c22*dy/R )
      !endif
    enddo    

  end subroutine mls2DDxSAThesis
!!---------------------End mls2DDxSAThesis---------------------!!



!!-----------------------findInvSymm3x3------------------------!!
  subroutine findInvSymm3x3(A,AInv,ADet)
  implicit none

    real(kind=C_K2),intent(in)::A(3,3)
    real(kind=C_K2),intent(out)::AInv(3,3),ADet

    AInv(1,1)=-A(2,3)**2 + A(2,2)*A(3,3) 
    AInv(1,2)=A(1,3)*A(2,3) - A(1,2)*A(3,3) 
    AInv(1,3)=-A(1,3)*A(2,2) + A(1,2)*A(2,3) 

    AInv(2,1)=AInv(1,2) 
    AInv(2,2)=-A(1,3)**2 + A(1,1)*A(3,3) 
    AInv(2,3)=A(1,2)*A(1,3) - A(1,1)*A(2,3) 
   
    AInv(3,1)=AInv(1,3)
    AInv(3,2)=AInv(2,3)
    AInv(3,3)=-A(1,2)**2 + A(1,1)*A(2,2)

    ADet = -A(1,3)*A(1,3)*A(2,2) + 2d0*A(1,2)*A(1,3)*A(2,3) &
      - A(1,1)*A(2,3)*A(2,3) - A(1,2)*A(1,2)*A(3,3) &
      + A(1,1)*A(2,2)*A(3,3)

    if(ADet.ne.0d0)then
      AInv = AInv/ADet
    else
      AInv=0d0
    endif

  end subroutine findInvSymm3x3
!!---------------------End findInvSymm3x3----------------------!!



  subroutine matMul_V13_M33_V31(V13,A,V31,res)
  implicit none

    real(kind=C_K2),intent(in)::V13(3),V31(3),A(3,3)
    real(kind=C_K2),intent(out)::res

    res=(A(1,1)*V13(1) + A(2,1)*V13(2) + A(3,1)*V13(3))*V31(1) &
      + (A(1,2)*V13(1) + A(2,2)*V13(2) + A(3,2)*V13(3))*V31(2) &
      + (A(1,3)*V13(1) + A(2,3)*V13(2) + A(3,3)*V13(3))*V31(3)

  end subroutine matMul_V13_M33_V31

end module meshFreeMod
!!-------------------------End meshFreeMod-------------------------!!



