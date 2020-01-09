subroutine matrixSet2(npoinl,npoint,nelem,conn,ivl,ivq,linkl,&
	linkq,invJ,totDe,u,v,gNAdv,gGxElev,gGyElev,por,porTurX,porTurY,&
  presK,cnst,stoneCnst,presOn)
implicit none
  
  !N advection matrix
  !G elevation gradient matrices x,y

  integer(kind=4),intent(in)::npoinl,npoint,nelem,conn(nelem,6)
  integer(kind=4),intent(in)::ivl(0:npoint),linkl(ivl(0)*npoint)
  integer(kind=4),intent(in)::ivq(0:npoint),linkq(ivq(0)*npoint)
  integer(kind=4)::i,j,k,i2,j2,k2,n(6),gRow,gCol,lRow,lCol
  integer(kind=4)::nlinkl(ivl(0)),nlinkq(ivq(0))  

  real(kind=8),intent(in)::invJ(nelem,5),totDe(npoinl),cnst(10)
  real(kind=8),intent(in)::u(npoinl),v(npoinl)
  real(kind=8),intent(in)::por(npoinl),stoneCnst(10)
  real(kind=8),intent(out)::gNAdv(ivq(0)*npoint)
  real(kind=8),intent(out)::gGxElev(ivl(0)*npoint)
  real(kind=8),intent(out)::gGyElev(ivl(0)*npoint)  
  real(kind=8),intent(out)::porTurX(ivq(0)*npoint)
  real(kind=8),intent(out)::porTurY(ivq(0)*npoint)
  real(kind=8),intent(out)::presK(ivq(0)*npoint)
  real(kind=8)::nAd(6,6),gxElev(6,3),gyElev(6,3)
  real(kind=8)::porCell,stTur,lporTurX(6,6),lporTurY(6,6)
  real(kind=8)::lpresK(6,6),tempr(3)

  logical,intent(in)::presOn

  gNAdv=0d0
  gGxElev=0d0
  gGyElev=0d0
  porTurX=0d0
  porTurY=0d0
  presK=0d0

  lpresK=0d0

  do i=1,nelem
    n=conn(i,:)

    do j=1,3
      tempr(j)=por(n(j))*totDe(n(j))
    enddo

    call nAdvMat(nAd,invJ(i,1),invJ(i,2),invJ(i,3),invJ(i,4),&
      u(n(1)),u(n(2)),u(n(3)),v(n(1)),v(n(2)),v(n(3)))

    call gFluxMat(gxElev,gyElev,invJ(i,1),invJ(i,2),invJ(i,3),&
      invJ(i,4),totDe(n(1)),totDe(n(2)),totDe(n(3)))    

    call porTurCal(lporTurX,u(n(1)),u(n(2)),u(n(3)))
    call porTurCal(lporTurY,v(n(1)),v(n(2)),v(n(3)))

    if(presOn) then
      call presCal(lpresK,tempr)
    endif

    porCell=(por(n(1))+por(n(2))+por(n(3)))/3.0d0
    stTur=stoneCnst(6)*(1-porCell)/(porCell**2)

    nAd=nAd*invJ(i,5)/porCell
    gxElev=gxElev*cnst(1)*invJ(i,5)*porCell
    gyElev=gyElev*cnst(1)*invJ(i,5)*porCell
    lporTurX=lporTurX*invJ(i,5)*stTur
    lporTurY=lporTurY*invJ(i,5)*stTur
    lpresK=lpresK*invJ(i,5)/cnst(8)

    !6x6
    do lRow=1,6
      gRow=n(lRow)
      k=(gRow-1)*ivq(0)
      nlinkq=linkq(k+1:k+ivq(0))
      do lCol=1,6
        gCol=n(lCol)
        do j=1,ivq(gRow)
          if(nlinkq(j).eq.gCol) goto 11
        enddo
        write(9,*)"[Err] node conn missing in N Advec at",gRow
        stop
        11 gNAdv(k+j)=gNAdv(k+j)+nAd(lRow,lCol)
        porTurX(k+j)=porTurX(k+j)+lporTurX(lRow,lCol)
        porTurY(k+j)=porTurY(k+j)+lporTurY(lRow,lCol)
        presK(k+j)=presK(k+j)+lpresK(lRow,lCol)
      enddo
    enddo

    !6x3
    do lRow=1,6
      gRow=n(lRow)
      k=(gRow-1)*ivl(0)
      nlinkl=linkl(k+1:k+ivl(0))
      do lCol=1,3
        gCol=n(lCol)
        do j=1,ivl(gRow)
          if(nlinkl(j).eq.gCol) goto 12
        enddo
        write(9,*)"[Err] node conn missing in G Elev Mat at",gRow
        stop
        12 gGxElev(k+j)=gGxElev(k+j)+gxElev(lRow,lCol)
        gGyElev(k+j)=gGyElev(k+j)+gyElev(lRow,lCol)
      enddo
    enddo

  enddo

end subroutine matrixSet2

subroutine nAdvMat(nAd,b11,b12,b21,b22,u1,u2,u3,v1,v2,v3)
implicit none
  
  real(kind=8),intent(out)::nAd(6,6)
  real(kind=8),intent(in)::b11,b12,b21,b22,u1,u2,u3,v1,v2,v3

  nAd(1,1)=(1d0/120d0)*(b11*(-8d0*u1 + u2 - u3) + b12*(-8d0*u1 &
    - u2 + u3) - 8d0*b21*v1 - 8d0*b22*v1 + b21*v2 - b22*v2 &
    - b21*v3 + b22*v3)

  nAd(1,2)=(1d0/360d0)*(b12*(u1 - u3) - b11*(5d0*u1 + 6d0*u2 &
    + u3) - 5d0*b21*v1 + b22*v1 - 6d0*b21*v2 - b21*v3 - b22*v3)

  nAd(1,3)=(1d0/360d0)*(b11*(u1 - u2) - b12*(5d0*u1 + u2 &
    + 6d0*u3) + b21*v1 - 5d0*b22*v1 - b21*v2 - b22*v2 &
    - 6d0*b22*v3)

  nAd(1,4)=(1d0/90d0)*(b12*(2d0*u2 + u3) + b11*(6d0*u1 + 2d0*u2 &
    + u3) + 6d0*b21*v1 + 2d0*b21*v2 + 2d0*b22*v2 + b21*v3 &
    + b22*v3)

  nAd(1,5)=(1d0/90d0)*(b11*(u1 - 2d0*(u2 + u3)) + b12*(u1 &
    - 2d0*(u2 + u3)) + (b21 + b22)*(v1 - 2d0*(v2 + v3)))

  nAd(1,6)=(1d0/90d0)*(b11*(u2 + 2d0*u3) + b12*(6d0*u1 + u2 &
    + 2d0*u3) + 6d0*b22*v1 + b21*v2 + b22*v2 + 2d0*b21*v3 &
    + 2d0*b22*v3)

  nAd(2,1)=(1d0/360d0)*(6d0*b12*(u1 + u2) + b11*(6d0*u1 &
    + 5d0*u2 + u3) + 6d0*b21*v1 + 6d0*b22*v1 + 5d0*b21*v2 &
    + 6d0*b22*v2 + b21*v3)

  nAd(2,2)=(1d0/120d0)*(-2d0*b12*u1 + 2d0*b12*u3 + b11*(-u1 &
    + 8d0*u2 + u3) - b21*v1 - 2d0*b22*v1 + 8d0*b21*v2 &
    + b21*v3 + 2d0*b22*v3)

  nAd(2,3)=(1d0/360d0)*(b11*(u1 - u2) - 6d0*b12*(u2 + u3) &
    + b21*v1 - b21*v2 - 6d0*b22*v2 - 6d0*b22*v3)

  nAd(2,4)=(1d0/90d0)*(-6d0*b12*u2 - b11*(2d0*u1 + 6d0*u2 &
    + u3) - 2d0*b21*v1 - 6d0*b21*v2 - 6d0*b22*v2 - b21*v3)

  nAd(2,5)=(1d0/90d0)*(6d0*b12*u2 - b11*(u1 + 2d0*u3) &
    - b21*v1 + 6d0*b22*v2 - 2d0*b21*v3)

  nAd(2,6)=(1d0/90d0)*(b11*(2d0*u1 - u2 + 2d0*u3) &
    + b21*(2d0*v1 - v2 + 2d0*v3))

  nAd(3,1)=(1d0/360d0)*(6d0*b11*(u1 + u3) + b12*(6d0*u1 + u2 &
    + 5d0*u3) + 6d0*b21*v1 + 6d0*b22*v1 + b22*v2 + 6d0*b21*v3 &
    + 5d0*b22*v3)

  nAd(3,2)=(1d0/360d0)*(b12*(u1 - u3) - 6d0*b11*(u2 + u3) &
    + b22*v1 - 6d0*b21*v2 - 6d0*b21*v3 - b22*v3)

  nAd(3,3)=(1d0/120d0)*(-2d0*b11*(u1 - u2) + b12*(-u1 + u2 &
    + 8d0*u3) - 2d0*b21*v1 - b22*v1 + 2d0*b21*v2 + b22*v2 &
    + 8d0*b22*v3)

  nAd(3,4)=(1d0/90d0)*(b12*(2d0*u1 + 2d0*u2 - u3) &
    + b22*(2d0*v1 + 2d0*v2 - v3))

  nAd(3,5)=(1d0/90d0)*((-b12)*(u1 + 2d0*u2) + 6d0*b11*u3 &
    - b22*v1 - 2d0*b22*v2 + 6d0*b21*v3)

  nAd(3,6)=(1d0/90d0)*(-6d0*b11*u3 - b12*(2d0*u1 + u2 &
    + 6d0*u3) - 2d0*b22*v1 - b22*v2 - 6d0*b21*v3 - 6d0*b22*v3)

  nAd(4,1)=(1d0/90d0)*((-(b11 + b12))*(6d0*u1 + 2d0*u2 &
    + u3) - (b21 + b22)*(6d0*v1 + 2d0*v2 + v3))

  nAd(4,2)=(1d0/90d0)*(b11*(2d0*u1 + 6d0*u2 + u3) &
    + b21*(2d0*v1 + 6d0*v2 + v3))

  nAd(4,3)=(1d0/90d0)*(b11*(u1 - u2) - b12*(u1 + 2d0*u2) &
    + b21*v1 - b22*v1 - b21*v2 - 2d0*b22*v2)

  nAd(4,4)=(-(2d0/45d0))*(b11*(u1 - u2) + b12*(4d0*u1 + 3d0*u2 &
    - u3) + b21*v1 + 4d0*b22*v1 - b21*v2 + 3d0*b22*v2 - b22*v3)

  nAd(4,5)=(2d0/45d0)*(b11*(2d0*u2 + u3) + b12*(u1 + 3d0*u2 &
    + 2d0*u3) + b22*v1 + 2d0*b21*v2 + 3d0*b22*v2 + b21*v3 &
    + 2d0*b22*v3)

  nAd(4,6)=(2d0/45d0)*((-b11)*(2d0*u1 + u3) + b12*(u1 + u2 &
    + u3) - 2d0*b21*v1 + b22*v1 + b22*v2 - b21*v3 + b22*v3)

  nAd(5,1)=(1d0/90d0)*(b12*(2d0*u2 + u3) + b11*(u2 + 2d0*u3) &
    + b21*v2 + 2d0*b22*v2 + 2d0*b21*v3 + b22*v3)

  nAd(5,2)=(1d0/90d0)*(b11*(u1 + 6d0*u2 + 2d0*u3) + b21*(v1 &
    + 6d0*v2 + 2d0*v3))

  nAd(5,3)=(1d0/90d0)*(b12*(u1 + 2d0*u2 + 6d0*u3) + b22*(v1 &
    + 2d0*v2 + 6d0*v3))

  nAd(5,4)=(-(2d0/45d0))*(b11*(u1 + u2 + u3) + b12*(2d0*u1 &
    + 3d0*u2 + u3) + b21*v1 + 2d0*b22*v1 + b21*v2 &
    + 3d0*b22*v2 + b21*v3 + b22*v3)

  nAd(5,5)=(-(2d0/45d0))*(b12*(u1 - 3d0*u2 - 4d0*u3) &
    + b11*(u1 - 4d0*u2 - 3d0*u3) + b21*v1 + b22*v1 &
    - 4d0*b21*v2 - 3d0*b22*v2 - 3d0*b21*v3 - 4d0*b22*v3)

  nAd(5,6)=(-(2d0/45d0))*(b12*(u1 + u2 + u3) + b11*(2d0*u1 &
    + u2 + 3d0*u3) + 2d0*b21*v1 + b22*v1 + b21*v2 + b22*v2 &
    + 3d0*b21*v3 + b22*v3)

  nAd(6,1)=(1d0/90d0)*((-(b11 + b12))*(6d0*u1 + u2 + 2d0*u3) &
    - (b21 + b22)*(6d0*v1 + v2 + 2d0*v3))

  nAd(6,2)=(1d0/90d0)*(b12*(u1 - u3) - b11*(u1 + 2d0*u3) &
    - b21*v1 + b22*v1 - 2d0*b21*v3 - b22*v3)

  nAd(6,3)=(1d0/90d0)*(b12*(2d0*u1 + u2 + 6d0*u3) + b22*(2d0*v1 &
    + v2 + 6d0*v3))

  nAd(6,4)=(2d0/45d0)*((-b12)*(2d0*u1 + u2) + b11*(u1 + u2 &
    + u3) + b21*v1 - 2d0*b22*v1 + b21*v2 - b22*v2 + b21*v3)

  nAd(6,5)=(2d0/45d0)*(b12*(u2 + 2d0*u3) + b11*(u1 + 2d0*u2 &
    + 3d0*u3) + b21*v1 + 2d0*b21*v2 + b22*v2 + 3d0*b21*v3 &
    + 2d0*b22*v3)

  nAd(6,6)=(-(2d0/45d0))*(b12*(u1 - u3) + b11*(4d0*u1 - u2 &
    + 3d0*u3) + 4d0*b21*v1 + b22*v1 - b21*v2 + 3d0*b21*v3 &
    - b22*v3)

end subroutine nAdvMat

subroutine gFluxMat(gxElev,gyElev,b11,b12,b21,b22,d1,d2,d3)
implicit none

  real(kind=8),intent(out)::gxElev(6,3),gyElev(6,3)
  real(kind=8),intent(in)::b11,b12,b21,b22,d1,d2,d3

  gxElev(1,1)=(-(1d0/120d0))*(b11 + b12)*(2d0*d1 - d2 - d3)

  gxElev(1,2)=(1d0/120d0)*b11*(2d0*d1 - d2 - d3)

  gxElev(1,3)=(1d0/120d0)*b12*(2d0*d1 - d2 - d3)

  gxElev(2,1)=(1d0/120d0)*(b11 + b12)*(d1 - 2d0*d2 + d3)

  gxElev(2,2)=(-(1d0/120d0))*b11*(d1 - 2d0*d2 + d3)

  gxElev(2,3)=(-(1d0/120d0))*b12*(d1 - 2d0*d2 + d3)

  gxElev(3,1)=(1d0/120d0)*(b11 + b12)*(d1 + d2 - 2d0*d3)

  gxElev(3,2)=(-(1d0/120d0))*b11*(d1 + d2 - 2d0*d3)

  gxElev(3,3)=(-(1d0/120d0))*b12*(d1 + d2 - 2d0*d3)

  gxElev(4,1)=(-(1d0/30d0))*(b11 + b12)*(2d0*d1 + 2d0*d2 + d3)

  gxElev(4,2)=(1d0/30d0)*b11*(2d0*d1 + 2d0*d2 + d3)

  gxElev(4,3)=(1d0/30d0)*b12*(2d0*d1 + 2d0*d2 + d3)

  gxElev(5,1)=(-(1d0/30d0))*(b11 + b12)*(d1 + 2d0*(d2 + d3))

  gxElev(5,2)=(1d0/30d0)*b11*(d1 + 2d0*(d2 + d3))

  gxElev(5,3)=(1d0/30d0)*b12*(d1 + 2d0*(d2 + d3))

  gxElev(6,1)=(-(1d0/30d0))*(b11 + b12)*(2d0*d1 + d2 + 2d0*d3)

  gxElev(6,2)=(1d0/30d0)*b11*(2d0*d1 + d2 + 2d0*d3)

  gxElev(6,3)=(1d0/30d0)*b12*(2d0*d1 + d2 + 2d0*d3)

  gyElev(1,1)=(-(1d0/120d0))*(b21 + b22)*(2d0*d1 - d2 - d3)

  gyElev(1,2)=(1d0/120d0)*b21*(2d0*d1 - d2 - d3)

  gyElev(1,3)=(1d0/120d0)*b22*(2d0*d1 - d2 - d3)

  gyElev(2,1)=(1d0/120d0)*(b21 + b22)*(d1 - 2d0*d2 + d3)

  gyElev(2,2)=(-(1d0/120d0))*b21*(d1 - 2d0*d2 + d3)

  gyElev(2,3)=(-(1d0/120d0))*b22*(d1 - 2d0*d2 + d3)

  gyElev(3,1)=(1d0/120d0)*(b21 + b22)*(d1 + d2 - 2d0*d3)

  gyElev(3,2)=(-(1d0/120d0))*b21*(d1 + d2 - 2d0*d3)

  gyElev(3,3)=(-(1d0/120d0))*b22*(d1 + d2 - 2d0*d3)

  gyElev(4,1)=(-(1d0/30d0))*(b21 + b22)*(2d0*d1 + 2d0*d2 + d3)

  gyElev(4,2)=(1d0/30d0)*b21*(2d0*d1 + 2d0*d2 + d3)

  gyElev(4,3)=(1d0/30d0)*b22*(2d0*d1 + 2d0*d2 + d3)

  gyElev(5,1)=(-(1d0/30d0))*(b21 + b22)*(d1 + 2d0*(d2 + d3))

  gyElev(5,2)=(1d0/30d0)*b21*(d1 + 2d0*(d2 + d3))

  gyElev(5,3)=(1d0/30d0)*b22*(d1 + 2d0*(d2 + d3))

  gyElev(6,1)=(-(1d0/30d0))*(b21 + b22)*(2d0*d1 + d2 + 2d0*d3)

  gyElev(6,2)=(1d0/30d0)*b21*(2d0*d1 + d2 + 2d0*d3)

  gyElev(6,3)=(1d0/30d0)*b22*(2d0*d1 + d2 + 2d0*d3)

end subroutine gFluxMat

! Modifed integration
! subroutine gFluxMat(gxElev,gyElev,b11,b12,b21,b22,d1,d2,d3)
! implicit none

!   real(kind=8),intent(out)::gxElev(6,3),gyElev(6,3)
!   real(kind=8),intent(in)::b11,b12,b21,b22,d1,d2,d3

!   gxElev(1,1)=(1d0/120d0)*(b12*(16*d1 + 3*d2 + d3) &
!     + b11*(16*d1 + d2 + 3*d3))
 
!   gxElev(1,2)=(1d0/120d0)*(2*b12*(d1 - d2) + b11*(2*d1 - d2 &
!     - d3))
   
!   gxElev(1,3)=(1d0/120d0)*(2*b11*(d1 - d3) + b12*(2*d1 &
!     - d2 - d3))
   
!   gxElev(2,1)=(1d0/120d0)*(b12*(-d1 + d3) + b11*(d1 &
!     - 2*d2 + d3))
   
!   gxElev(2,2)=(1d0/120d0)*(2*b12*(d1 - d3) - b11*(d1 &
!     + 16*d2 + 3*d3))
   
!   gxElev(2,3)=(1d0/120d0)*(b12*(-d1 + d3) &
!     + 2*b11*(-d2 + d3))
   
!   gxElev(3,1)=(1d0/120d0)*(b11*(-d1 + d2) + b12*(d1 &
!     + d2 - 2*d3))
   
!   gxElev(3,2)=(1d0/120d0)*(b11*(-d1 + d2) &
!     + 2*b12*(d2 - d3))
   
!   gxElev(3,2)=(1d0/120d0)*(2*b11*(d1 - d2) - b12*(d1 &
!     + 3*d2 + 16*d3))
   
!   gxElev(4,1)=(1d0/30d0)*(b12*(4*d1 + 2*d2 - d3) &
!     - b11*(2*d1 + 2*d2 + d3))
   
!   gxElev(4,2)=(1d0/30d0)*(b12*(4*d1 + 6*d2) + b11*(2*d1 &
!     + 2*d2 + d3))
   
!   gxElev(4,3)=(1d0/30d0)*b12*(2*d1 + 2*d2 + d3)
   
!   gxElev(5,1)=(-(1d0/30d0))*(b11 + b12)*(d1 + 2*(d2 + d3))
   
!   gxElev(5,2)=(1d0/30d0)*(-2*b12*(3*d2 + 2*d3) &
!     + b11*(d1 - 2*(2*d2 + d3)))
   
!   gxElev(5,3)=(1d0/30d0)*(-2*b11*(2*d2 + 3*d3) &
!     + b12*(d1 - 2*(d2 + 2*d3)))
   
!   gxElev(6,1)=(1d0/30d0)*(b11*(4*d1 - d2 + 2*d3) &
!     - b12*(2*d1 + d2 + 2*d3))
   
!   gxElev(6,2)=(1d0/30d0)*b11*(2*d1 + d2 + 2*d3)
   
!   gxElev(6,3)=(1d0/30d0)*(b12*(2*d1 + d2 + 2*d3) &
!     + b11*(4*d1 + 6*d3))


!   gyElev(1,1)=(1d0/120d0)*(b22*(16*d1 + 3*d2 + d3) &
!     + b21*(16*d1 + d2 + 3*d3))
 
!   gyElev(1,2)=(1d0/120d0)*(2*b22*(d1 - d2) + b21*(2*d1 - d2 &
!     - d3))
   
!   gyElev(1,3)=(1d0/120d0)*(2*b21*(d1 - d3) + b22*(2*d1 &
!     - d2 - d3))
   
!   gyElev(2,1)=(1d0/120d0)*(b22*(-d1 + d3) + b21*(d1 &
!     - 2*d2 + d3))
   
!   gyElev(2,2)=(1d0/120d0)*(2*b22*(d1 - d3) - b21*(d1 &
!     + 16*d2 + 3*d3))
   
!   gyElev(2,3)=(1d0/120d0)*(b22*(-d1 + d3) &
!     + 2*b21*(-d2 + d3))
   
!   gyElev(3,1)=(1d0/120d0)*(b21*(-d1 + d2) + b22*(d1 &
!     + d2 - 2*d3))
   
!   gyElev(3,2)=(1d0/120d0)*(b21*(-d1 + d2) &
!     + 2*b22*(d2 - d3))
   
!   gyElev(3,2)=(1d0/120d0)*(2*b21*(d1 - d2) - b22*(d1 &
!     + 3*d2 + 16*d3))
   
!   gyElev(4,1)=(1d0/30d0)*(b22*(4*d1 + 2*d2 - d3) &
!     - b21*(2*d1 + 2*d2 + d3))
   
!   gyElev(4,2)=(1d0/30d0)*(b22*(4*d1 + 6*d2) + b21*(2*d1 &
!     + 2*d2 + d3))
   
!   gyElev(4,3)=(1d0/30d0)*b22*(2*d1 + 2*d2 + d3)
   
!   gyElev(5,1)=(-(1d0/30d0))*(b21 + b22)*(d1 + 2*(d2 + d3))
   
!   gyElev(5,2)=(1d0/30d0)*(-2*b22*(3*d2 + 2*d3) &
!     + b21*(d1 - 2*(2*d2 + d3)))
   
!   gyElev(5,3)=(1d0/30d0)*(-2*b21*(2*d2 + 3*d3) &
!     + b22*(d1 - 2*(d2 + 2*d3)))
   
!   gyElev(6,1)=(1d0/30d0)*(b21*(4*d1 - d2 + 2*d3) &
!     - b22*(2*d1 + d2 + 2*d3))
   
!   gyElev(6,2)=(1d0/30d0)*b21*(2*d1 + d2 + 2*d3)
   
!   gyElev(6,3)=(1d0/30d0)*(b22*(2*d1 + d2 + 2*d3) &
!     + b21*(4*d1 + 6*d3))

! end subroutine gFluxMat

subroutine porTurCal(lporTur,u1,u2,u3)
implicit none

  real(kind=8),intent(out)::lporTur(6,6)
  real(kind=8),intent(in)::u1,u2,u3

  lporTur(1,1)=(1d0/420d0)*(5*u1 + u2 + u3)
  lporTur(1,2)=(-4*u1 - 4*u2 + u3)/2520d0
  lporTur(1,3)=(-4*u1 + u2 - 4*u3)/2520d0
  lporTur(1,4)=(1d0/630d0)*(3*u1 - 2*u2 - u3)
  lporTur(1,5)=(1d0/630d0)*(-u1 - 3*(u2 + u3))
  lporTur(1,6)=(1d0/630d0)*(3*u1 - u2 - 2*u3)


  lporTur(2,1)=(-4*u1 - 4*u2 + u3)/2520d0
  lporTur(2,2)=(1d0/420d0)*(u1 + 5*u2 + u3)
  lporTur(2,3)=(u1 - 4*(u2 + u3))/2520d0
  lporTur(2,4)=(1d0/630d0)*(-2*u1 + 3*u2 - u3)
  lporTur(2,5)=(1d0/630d0)*(-u1 + 3*u2 - 2*u3)
  lporTur(2,6)=(1d0/630d0)*(-3*u1 - u2 - 3*u3)


  lporTur(3,1)=(-4*u1 + u2 - 4*u3)/2520d0
  lporTur(3,2)=(u1 - 4*(u2 + u3))/2520d0
  lporTur(3,3)=(1d0/420d0)*(u1 + u2 + 5*u3)
  lporTur(3,4)=(1d0/630d0)*(-3*u1 - 3*u2 - u3)
  lporTur(3,5)=(1d0/630d0)*(-u1 - 2*u2 + 3*u3)
  lporTur(3,6)=(1d0/630d0)*(-2*u1 - u2 + 3*u3)


  lporTur(4,1)=(1d0/630d0)*(3*u1 - 2*u2 - u3)
  lporTur(4,2)=(1d0/630d0)*(-2*u1 + 3*u2 - u3)
  lporTur(4,3)=(1d0/630d0)*(-3*u1 - 3*u2 - u3)
  lporTur(4,4)=(4d0/315d0)*(3*u1 + 3*u2 + u3)
  lporTur(4,5)=(2d0/315d0)*(2*u1 + 3*u2 + 2*u3)
  lporTur(4,6)=(2d0/315d0)*(3*u1 + 2*(u2 + u3))

   
  lporTur(5,1)=(1d0/630d0)*(-u1 - 3*(u2 + u3))
  lporTur(5,2)=(1d0/630d0)*(-u1 + 3*u2 - 2*u3)
  lporTur(5,3)=(1d0/630d0)*(-u1 - 2*u2 + 3*u3)
  lporTur(5,4)=(2d0/315d0)*(2*u1 + 3*u2 + 2*u3)
  lporTur(5,5)=(4d0/315d0)*(u1 + 3*(u2 + u3))
  lporTur(5,6)=(2d0/315d0)*(2*u1 + 2*u2 + 3*u3)


  lporTur(6,1)=(1d0/630d0)*(3*u1 - u2 - 2*u3)
  lporTur(6,2)=(1d0/630d0)*(-3*u1 - u2 - 3*u3)
  lporTur(6,3)=(1d0/630d0)*(-2*u1 - u2 + 3*u3)
  lporTur(6,4)=(2d0/315d0)*(3*u1 + 2*(u2 + u3))
  lporTur(6,5)=(2d0/315d0)*(2*u1 + 2*u2 + 3*u3)
  lporTur(6,6)=(4d0/315d0)*(3*u1 + u2 + 3*u3)

end subroutine porTurCal

subroutine presCal(lpresK,q)
implicit none

  real(kind=8),intent(out)::lpresK(6,6)
  real(kind=8),intent(in)::q(3)

  lpresK(1,1)=(1d0/420d0)*(5*q(1) + q(2) + q(3))
  lpresK(1,2)=(-4*q(1) - 4*q(2) + q(3))/2520d0
  lpresK(1,3)=(-4*q(1) + q(2) - 4*q(3))/2520d0
  lpresK(1,4)=(1d0/630d0)*(3*q(1) - 2*q(2) - q(3))
  lpresK(1,5)=(1d0/630d0)*(-q(1) - 3*(q(2) + q(3)))
  lpresK(1,6)=(1d0/630d0)*(3*q(1) - q(2) - 2*q(3))

  lpresK(2,1)=(-4*q(1) - 4*q(2) + q(3))/2520d0
  lpresK(2,2)=(1d0/420d0)*(q(1) + 5*q(2) + q(3))
  lpresK(2,3)=(q(1) - 4*(q(2) + q(3)))/2520d0
  lpresK(2,4)=(1d0/630d0)*(-2*q(1) + 3*q(2) - q(3))
  lpresK(2,5)=(1d0/630d0)*(-q(1) + 3*q(2) - 2*q(3))
  lpresK(2,6)=(1d0/630d0)*(-3*q(1) - q(2) - 3*q(3))

  lpresK(3,1)=(-4*q(1) + q(2) - 4*q(3))/2520d0
  lpresK(3,2)=(q(1) - 4*(q(2) + q(3)))/2520d0
  lpresK(3,3)=(1d0/420d0)*(q(1) + q(2) + 5*q(3))
  lpresK(3,4)=(1d0/630d0)*(-3*q(1) - 3*q(2) - q(3))
  lpresK(3,5)=(1d0/630d0)*(-q(1) - 2*q(2) + 3*q(3))
  lpresK(3,6)=(1d0/630d0)*(-2*q(1) - q(2) + 3*q(3))

  lpresK(4,1)=(1d0/630d0)*(3*q(1) - 2*q(2) - q(3))
  lpresK(4,2)=(1d0/630d0)*(-2*q(1) + 3*q(2) - q(3))
  lpresK(4,3)=(1d0/630d0)*(-3*q(1) - 3*q(2) - q(3))
  lpresK(4,4)=(4d0/315d0)*(3*q(1) + 3*q(2) + q(3))
  lpresK(4,5)=(2d0/315d0)*(2*q(1) + 3*q(2) + 2*q(3))
  lpresK(4,6)=(2d0/315d0)*(3*q(1) + 2*(q(2) + q(3)))

  lpresK(5,1)=(1d0/630d0)*(-q(1) - 3*(q(2) + q(3)))
  lpresK(5,2)=(1d0/630d0)*(-q(1) + 3*q(2) - 2*q(3))
  lpresK(5,3)=(1d0/630d0)*(-q(1) - 2*q(2) + 3*q(3))
  lpresK(5,4)=(2d0/315d0)*(2*q(1) + 3*q(2) + 2*q(3))
  lpresK(5,5)=(4d0/315d0)*(q(1) + 3*(q(2) + q(3)))
  lpresK(5,6)=(2d0/315d0)*(2*q(1) + 2*q(2) + 3*q(3))

  lpresK(6,1)=(1d0/630d0)*(3*q(1) - q(2) - 2*q(3))
  lpresK(6,2)=(1d0/630d0)*(-3*q(1) - q(2) - 3*q(3))
  lpresK(6,3)=(1d0/630d0)*(-2*q(1) - q(2) + 3*q(3))
  lpresK(6,4)=(2d0/315d0)*(3*q(1) + 2*(q(2) + q(3)))
  lpresK(6,5)=(2d0/315d0)*(2*q(1) + 2*q(2) + 3*q(3))
  lpresK(6,6)=(4d0/315d0)*(3*q(1) + q(2) + 3*q(3))

end subroutine presCal