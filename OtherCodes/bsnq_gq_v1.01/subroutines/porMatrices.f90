subroutine porMatrices1(npoinl,npoint,nelem,conn,ivl,ivq,linkl,&
  linkq,invJ,depth,por,porLam,gCxFlux,gCyFlux,gAx,gAy,&
  cnst,stoneCnst)
implicit none
  
  integer(kind=4),intent(in)::npoinl,npoint,nelem,conn(nelem,6)
  integer(kind=4),intent(in)::ivl(0:npoint),linkl(ivl(0)*npoint)
  integer(kind=4),intent(in)::ivq(0:npoint),linkq(ivq(0)*npoint)
  integer(kind=4)::i,j,k,i2,j2,k2,n(6),gRow,gCol,lRow,lCol
  integer(kind=4)::nlinkl(ivl(0)),nlinkq(ivq(0))

  real(kind=8),intent(in)::invJ(nelem,5),depth(npoint),cnst(10)
  real(kind=8),intent(in)::por(npoinl),stoneCnst(10)
  real(kind=8),intent(out)::porLam(ivq(0)*npoint)  
  real(kind=8),intent(out)::gCxFlux(ivq(0)*npoinl)
  real(kind=8),intent(out)::gCyFlux(ivq(0)*npoinl)
  real(kind=8),intent(out)::gAx(ivl(0)*npoint),gAy(ivl(0)*npoint)
  real(kind=8)::cxFlux(3,6),cyFlux(3,6),locM1(6,6)
  real(kind=8)::ax(6,3),ay(6,3),porCell

  porLam=0d0
  gCxFlux=0d0
  gCyFlux=0d0
  gAx=0d0
  gAy=0d0

  locM1=0d0
  locM1(1,:)=(/ 1d0/60d0,-1d0/360d0,-1d0/360d0,0d0,-1d0/90d0,0d0 /)
  locM1(2,:)=(/ -1d0/360d0,1d0/60d0,-1d0/360d0,0d0,0d0,-1d0/90d0 /)
  locM1(3,:)=(/ -1d0/360d0,-1d0/360d0,1d0/60d0,-1d0/90d0,0d0,0d0 /)
  locM1(4,:)=(/ 0d0,0d0,-1d0/90d0,4d0/45d0,2d0/45d0,2d0/45d0 /)
  locM1(5,:)=(/ -1d0/90d0,0d0,0d0,2d0/45d0,4d0/45d0,2d0/45d0 /)
  locM1(6,:)=(/ 0d0,-1d0/90d0,0d0,2d0/45d0,2d0/45d0,4d0/45d0 /)

  do i=1,nelem    !go to each ele
    n=conn(i,:)  !get node conn
    porCell=(por(n(1))+por(n(2))+por(n(3)))/3d0    
    porCell=((1d0-porCell)**3d0)/porCell    
    do lRow=1,6   !iter thu local row
      gRow=n(lRow)     !global row
      k=(gRow-1)*ivq(0) !global row pos in 1d global stiff
      nlinkq=linkq(k+1:k+ivq(0))  !link for grow
      do lCol=1,6         !iter thru local col
        gCol=n(lCol)     !global col
        do j=1,ivq(gRow)  !search for gCol in link of gRow
          if(nlinkq(j).eq.gCol) goto 11
        enddo
        write(9,*)"[Err] node conn missing in Mass1 at",gRow
        stop
        11 porLam(k+j)=porLam(k+j)+(locM1(lRow,lCol)*invJ(i,5) &
            *stoneCnst(5)*porCell)                
      enddo
    enddo
  enddo

  do i=1,nelem
    n=conn(i,:)    
    
    call cFluxMatPor(cxFlux,cyFlux,invJ(i,1),invJ(i,2),&
      invJ(i,3),invJ(i,4))

    call remainMatPor(ax,ay,invJ(i,1),invJ(i,2),invJ(i,3),&
      invJ(i,4),depth(n(1)),depth(n(2)),depth(n(3)))    

    
    porCell=(por(n(1))+por(n(2))+por(n(3)))/3.0d0
    
    cxFlux=-invJ(i,5)*cxFlux/porCell
    cyFlux=-invJ(i,5)*cyFlux/porCell
    ax=invJ(i,5)*cnst(6)*ax*porCell
    ay=invJ(i,5)*cnst(6)*ay*porCell   

    !3x6    
    do lRow=1,3
      gRow=n(lRow)
      k=(gRow-1)*ivq(0)
      nlinkq=linkq(k+1:k+ivq(0))
      do lCol=1,6
        gCol=n(lCol)
        do j=1,ivq(gRow)
          if(nlinkq(j).eq.gCol) goto 12
        enddo
        write(9,*)"[Err] node conn missing in C Flux at",gRow
        stop
        12 gCxFlux(k+j)=gCxFlux(k+j)+cxFlux(lRow,lCol)
        gCyFlux(k+j)=gCyFlux(k+j)+cyFlux(lRow,lCol)
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
          if(nlinkl(j).eq.gCol) goto 13
        enddo
        write(9,*)"[Err] node conn missing in A Mat at",gRow
        stop
        13 gAx(k+j)=gAx(k+j)+ax(lRow,lCol)
        gAy(k+j)=gAy(k+j)+ay(lRow,lCol)
      enddo
    enddo

  enddo

end subroutine porMatrices1

subroutine cFluxMatPor(cxFlux,cyFlux,b11,b12,b21,b22)
implicit none

  real(kind=8),intent(out)::cxFlux(3,6),cyFlux(3,6)
  real(kind=8),intent(in)::b11,b12,b21,b22

  ! !As per paper dN3idx N6j
  ! cxFlux=0d0
  ! cxFlux(1,4)=(1d0/6d0)*(-b11 - b12)
  ! cxFlux(1,5)=(1d0/6d0)*(-b11 - b12)
  ! cxFlux(1,6)=(1d0/6d0)*(-b11 - b12)
  ! cxFlux(2,4)=b11/6d0
  ! cxFlux(2,5)=b11/6d0
  ! cxFlux(2,6)=b11/6d0
  ! cxFlux(3,4)=b12/6d0
  ! cxFlux(3,5)=b12/6d0
  ! cxFlux(3,6)=b12/6d0

  ! !As per paper dN3idy N6j
  ! cyFlux=0d0
  ! cyFlux(1,4)=(1d0/6d0)*(-b21 - b22)
  ! cyFlux(1,5)=(1d0/6d0)*(-b21 - b22)
  ! cyFlux(1,6)=(1d0/6d0)*(-b21 - b22)
  ! cyFlux(2,4)=b21/6d0
  ! cyFlux(2,5)=b21/6d0
  ! cyFlux(2,6)=b21/6d0
  ! cyFlux(3,4)=b22/6d0
  ! cyFlux(3,5)=b22/6d0
  ! cyFlux(3,6)=b22/6d0

  ! As per Shagun N3i dN6jdx
  cxFlux(1,1)=(1d0/6d0)*(-b11 - b12)
  cxFlux(1,2)=0d0
  cxFlux(1,3)=0d0
  cxFlux(1,4)=(b11 - b12)/6d0 
  cxFlux(1,5)=(b11 + b12)/6d0
  cxFlux(1,6)=(1d0/6d0)*(-b11 + b12)

  cxFlux(2,1)=0d0 
  cxFlux(2,2)=b11/6d0 
  cxFlux(2,3)=0d0
  cxFlux(2,4)=(1d0/6d0)*(-b11 - 2d0*b12)
  cxFlux(2,5)=(1d0/6d0)*(b11 + 2d0*b12) 
  cxFlux(2,6)=-(b11/6d0) 

  cxFlux(3,1)=0d0 
  cxFlux(3,2)=0d0 
  cxFlux(3,3)=b12/6d0 
  cxFlux(3,4)=-(b12/6d0)
  cxFlux(3,5)=(1d0/6d0)*(2d0*b11 + b12)
  cxFlux(3,6)=(1d0/6d0)*(-2d0*b11 - b12)

  ! As per Shagun N3i dN6jdy
  cyFlux(1,1)=(1d0/6d0)*(-b21 - b22)
  cyFlux(1,2)=0d0
  cyFlux(1,3)=0d0
  cyFlux(1,4)=(b21 - b22)/6d0 
  cyFlux(1,5)=(b21 + b22)/6d0
  cyFlux(1,6)=(1d0/6d0)*(-b21 + b22)

  cyFlux(2,1)=0d0
  cyFlux(2,2)=b21/6d0
  cyFlux(2,3)=0d0
  cyFlux(2,4)=(1d0/6d0)*(-b21 - 2d0*b22) 
  cyFlux(2,5)=(1d0/6d0)*(b21 + 2d0*b22)
  cyFlux(2,6)=-(b21/6d0)

  cyFlux(3,1)=0d0
  cyFlux(3,2)=0d0
  cyFlux(3,3)=b22/6d0 
  cyFlux(3,4)=-(b22/6d0)
  cyFlux(3,5)=(1d0/6d0)*(2d0*b21 + b22)
  cyFlux(3,6)=(1d0/6d0)*(-2d0*b21 - b22)

end subroutine cFluxMatPor

subroutine remainMatPor(ax,ay,b11,b12,b21,b22,h1,h2,h3)
implicit none

  real(kind=8),intent(out)::ax(6,3),ay(6,3)
  real(kind=8),intent(in)::b11,b12,b21,b22,h1,h2,h3 

  ax(1,1)=(1d0/180d0)*((-b12)*(21d0*h1**2d0 + 6d0*h1*h2 &
    + h2**2d0 + h2*h3 + h3**2d0) - b11*(21d0*h1**2d0 &
    + h2**2d0 + 6d0*h1*h3 + h2*h3 + h3**2d0))

  ax(1,2)=(1d0/180d0)*(3d0*b12*(-h1**2d0 + h2**2d0) &
    + b11*(-3d0*h1**2d0 + h2**2d0 + h2*h3 + h3**2d0))

  ax(1,3)=(1d0/180d0)*(-3d0*b11*(h1**2d0 - h3**2d0) &
    + b12*(-3d0*h1**2d0 + h2**2d0 + h2*h3 + h3**2d0))

  ax(2,1)=(1d0/180d0)*(b12*(2d0*h1**2d0 - h1*h3 - h3**2d0) &
    - b11*(h1**2d0 - 3d0*h2**2d0 + h1*h3 + h3**2d0))

  ax(2,2)=(1d0/180d0)*(6d0*b12*h2*(-h1 + h3) + b11*(h1**2d0 &
    + 21d0*h2**2d0 + h1*h3 + 6d0*h2*h3 + h3**2d0))

  ax(2,3)=(1d0/180d0)*(b12*(h1**2d0 + h1*h3 - 2d0*h3**2d0) &
    + 3d0*b11*(h2**2d0 - h3**2d0))

  ax(3,1)=(1d0/180d0)*(b11*(2d0*h1**2d0 - h1*h2 - h2**2d0) &
    - b12*(h1**2d0 + h1*h2 + h2**2d0 - 3d0*h3**2d0))

  ax(3,2)=(1d0/180d0)*(b11*(h1**2d0 + h1*h2 - 2d0*h2**2d0) &
    + 3d0*b12*(-h2**2d0 + h3**2d0))

  ax(3,3)=(1d0/180d0)*(6d0*b11*(-h1 + h2)*h3 + b12*(h1**2d0 &
    + h1*h2 + h2**2d0 + 6d0*h2*h3 + 21d0*h3**2d0))

  ax(4,1)=(1d0/90d0)*(b12*(-9d0*h1**2d0 - 8d0*h1*h2 &
    - 3d0*h2**2d0 + 2d0*h1*h3 + 2d0*h2*h3 + h3**2d0) &
    + b11*(3d0*h1**2d0 + 3d0*h2**2d0 + 2d0*h2*h3 + h3**2d0 &
    + 2d0*h1*(2d0*h2 + h3)))

  ax(4,2)=(1d0/90d0)*(-6d0*b12*(h1**2d0 + 2d0*h1*h2 &
    + 2d0*h2**2d0) - b11*(3d0*h1**2d0 + 3d0*h2**2d0 &
    + 2d0*h2*h3 + h3**2d0 + 2d0*h1*(2d0*h2 + h3)))

  ax(4,3)=(-(1d0/90d0))*b12*(3d0*h1**2d0 + 3d0*h2**2d0 &
    + 2d0*h2*h3 + h3**2d0 + 2d0*h1*(2d0*h2 + h3))

  ax(5,1)=(1d0/90d0)*(b11 + b12)*(h1**2d0 + 3d0*h2**2d0 &
    + 4d0*h2*h3 + 3d0*h3**2d0 + 2d0*h1*(h2 + h3))

  ax(5,2)=(1d0/90d0)*(6d0*b12*(2d0*h2**2d0 + 2d0*h2*h3 &
    + h3**2d0) - b11*(h1**2d0 - 9d0*h2**2d0 - 8d0*h2*h3 &
    - 3d0*h3**2d0 + 2d0*h1*(h2 + h3)))

  ax(5,3)=(1d0/90d0)*(6d0*b11*(h2**2d0 + 2d0*h2*h3 &
    + 2d0*h3**2d0) - b12*(h1**2d0 - 3d0*h2**2d0 - 8d0*h2*h3 &
    - 9d0*h3**2d0 + 2d0*h1*(h2 + h3)))

  ax(6,1)=(1d0/90d0)*(b11*(-9d0*h1**2d0 + h2**2d0 + 2d0*h1*(h2 &
    - 4d0*h3) + 2d0*h2*h3 - 3d0*h3**2d0) + b12*(3d0*h1**2d0 &
    + h2**2d0 + 2d0*h2*h3 + 3d0*h3**2d0 + 2d0*h1*(h2 + 2d0*h3)))

  ax(6,2)=(-(1d0/90d0))*b11*(3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 + 2d0*h1*(h2 + 2d0*h3))

  ax(6,3)=(1d0/90d0)*(-6d0*b11*(h1**2d0 + 2d0*h1*h3 &
    + 2d0*h3**2d0) - b12*(3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 + 2d0*h1*(h2 + 2d0*h3))) 


  ay(1,1)=(1d0/180d0)*((-b22)*(21d0*h1**2d0 + 6d0*h1*h2 &
    + h2**2d0 + h2*h3 + h3**2d0) - b21*(21d0*h1**2d0 &
    + h2**2d0 + 6d0*h1*h3 + h2*h3 + h3**2d0))

  ay(1,2)=(1d0/180d0)*(3d0*b22*(-h1**2d0 + h2**2d0) &
    + b21*(-3d0*h1**2d0 + h2**2d0 + h2*h3 + h3**2d0))

  ay(1,3)=(1d0/180d0)*(-3d0*b21*(h1**2d0 - h3**2d0) &
    + b22*(-3d0*h1**2d0 + h2**2d0 + h2*h3 + h3**2d0))

  ay(2,1)=(1d0/180d0)*(b22*(2d0*h1**2d0 - h1*h3 - h3**2d0) &
    - b21*(h1**2d0 - 3d0*h2**2d0 + h1*h3 + h3**2d0))

  ay(2,2)=(1d0/180d0)*(6d0*b22*h2*(-h1 + h3) + b21*(h1**2d0 &
    + 21d0*h2**2d0 + h1*h3 + 6d0*h2*h3 + h3**2d0))

  ay(2,3)=(1d0/180d0)*(b22*(h1**2d0 + h1*h3 - 2d0*h3**2d0) &
    + 3d0*b21*(h2**2d0 - h3**2d0))

  ay(3,1)=(1d0/180d0)*(b21*(2d0*h1**2d0 - h1*h2 - h2**2d0) &
    - b22*(h1**2d0 + h1*h2 + h2**2d0 - 3d0*h3**2d0))

  ay(3,2)=(1d0/180d0)*(b21*(h1**2d0 + h1*h2 - 2d0*h2**2d0) &
    + 3d0*b22*(-h2**2d0 + h3**2d0))

  ay(3,3)=(1d0/180d0)*(6d0*b21*(-h1 + h2)*h3 + b22*(h1**2d0 &
    + h1*h2 + h2**2d0 + 6d0*h2*h3 + 21d0*h3**2d0))

  ay(4,1)=(1d0/90d0)*(b22*(-9d0*h1**2d0 - 8d0*h1*h2 &
    - 3d0*h2**2d0 + 2d0*h1*h3 + 2d0*h2*h3 + h3**2d0) &
    + b21*(3d0*h1**2d0 + 3d0*h2**2d0 + 2d0*h2*h3 + h3**2d0 &
    + 2d0*h1*(2d0*h2 + h3)))

  ay(4,2)=(1d0/90d0)*(-6d0*b22*(h1**2d0 + 2d0*h1*h2 &
    + 2d0*h2**2d0) - b21*(3d0*h1**2d0 + 3d0*h2**2d0 &
    + 2d0*h2*h3 + h3**2d0 + 2d0*h1*(2d0*h2 + h3)))

  ay(4,3)=(-(1d0/90d0))*b22*(3d0*h1**2d0 + 3d0*h2**2d0 &
    + 2d0*h2*h3 + h3**2d0 + 2d0*h1*(2d0*h2 + h3))

  ay(5,1)=(1d0/90d0)*(b21 + b22)*(h1**2d0 + 3d0*h2**2d0 &
    + 4d0*h2*h3 + 3d0*h3**2d0 + 2d0*h1*(h2 + h3))

  ay(5,2)=(1d0/90d0)*(6d0*b22*(2d0*h2**2d0 + 2d0*h2*h3 &
    + h3**2d0) - b21*(h1**2d0 - 9d0*h2**2d0 - 8d0*h2*h3 &
    - 3d0*h3**2d0 + 2d0*h1*(h2 + h3)))

  ay(5,3)=(1d0/90d0)*(6d0*b21*(h2**2d0 + 2d0*h2*h3 &
    + 2d0*h3**2d0) - b22*(h1**2d0 - 3d0*h2**2d0 - 8d0*h2*h3 &
    - 9d0*h3**2d0 + 2d0*h1*(h2 + h3)))

  ay(6,1)=(1d0/90d0)*(b21*(-9d0*h1**2d0 + h2**2d0 &
    + 2d0*h1*(h2 - 4d0*h3) + 2d0*h2*h3 - 3d0*h3**2d0) &
    + b22*(3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 + 3d0*h3**2d0 &
    + 2d0*h1*(h2 + 2d0*h3)))

  ay(6,2)=(-(1d0/90d0))*b21*(3d0*h1**2d0 + h2**2d0 &
    + 2d0*h2*h3 + 3d0*h3**2d0 + 2d0*h1*(h2 + 2d0*h3))

  ay(6,3)=(1d0/90d0)*(-6d0*b21*(h1**2d0 + 2d0*h1*h3 &
    + 2d0*h3**2d0) - b22*(3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 + 2d0*h1*(h2 + 2d0*h3))) 


end subroutine remainMatPor
