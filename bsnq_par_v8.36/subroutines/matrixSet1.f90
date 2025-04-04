subroutine massMatrices(npoint,npoinl,nelem,conn,ivl,ivq,linkl,&
  linkq,invJ,mass1,mass2,porLam,por,stoneCnst,touMass)
implicit none
  
  integer(kind=4),intent(in)::npoint,npoinl,nelem
  integer(kind=4),intent(in)::ivl(0:npoint),ivq(0:npoint)
  integer(kind=4),intent(in)::linkl(ivl(0)*npoint)
  integer(kind=4),intent(in)::linkq(ivq(0)*npoint),conn(nelem,6)
  integer(kind=4)::i,j,k,i2,j2,k2,nq(6),nl(3),gRow,gCol,lRow,lCol
  integer(kind=4)::nlinkl(ivl(0)),nlinkq(ivq(0))

  real(kind=8),intent(in)::invJ(nelem,5),stoneCnst(10),por(npoinl)
  real(kind=8),intent(out)::porLam(ivq(0)*npoint)  
  real(kind=8),intent(out)::touMass(ivl(0)*npoint)
  real(kind=8),intent(out)::mass1(ivq(0)*npoint)
  real(kind=8),intent(out)::mass2(ivl(0)*npoinl)
  real(kind=8)::locM1(6,6),locM2(3,3),porCell
  real(kind=8)::loctouMass(6,3)
  
  locM1=0d0
  locM1(1,:)=(/ 1d0/60d0,-1d0/360d0,-1d0/360d0,0d0,-1d0/90d0,0d0 /)
  locM1(2,:)=(/ -1d0/360d0,1d0/60d0,-1d0/360d0,0d0,0d0,-1d0/90d0 /)
  locM1(3,:)=(/ -1d0/360d0,-1d0/360d0,1d0/60d0,-1d0/90d0,0d0,0d0 /)
  locM1(4,:)=(/ 0d0,0d0,-1d0/90d0,4d0/45d0,2d0/45d0,2d0/45d0 /)
  locM1(5,:)=(/ -1d0/90d0,0d0,0d0,2d0/45d0,4d0/45d0,2d0/45d0 /)
  locM1(6,:)=(/ 0d0,-1d0/90d0,0d0,2d0/45d0,2d0/45d0,4d0/45d0 /)

  locM2=0d0
  locM2(1,:)=(/ 1d0/12d0,1d0/24d0,1d0/24d0 /)
  locM2(2,:)=(/ 1d0/24d0,1d0/12d0,1d0/24d0 /)
  locM2(3,:)=(/ 1d0/24d0,1d0/24d0,1d0/12d0 /) 

  loctouMass=0d0
  loctouMass(1,:)=(/ 1d0/60d0,-1d0/120d0,-1d0/120d0 /)
  loctouMass(2,:)=(/ -1d0/120d0,1d0/60d0,-1d0/120d0 /)
  loctouMass(3,:)=(/ -1d0/120d0,-1d0/120d0,1d0/60d0 /)
  loctouMass(4,:)=(/ 1d0/15d0,1d0/15d0,1d0/30d0 /)
  loctouMass(5,:)=(/ 1d0/30d0,1d0/15d0,1d0/15d0 /)
  loctouMass(6,:)=(/ 1d0/15d0,1d0/30d0,1d0/15d0 /)


  ! do i=1,6
  !   write(9,*)locM1(i,:)
  ! enddo
  ! write(9,*)
  ! do i=1,3
  !   write(9,*)locM2(i,:)
  ! enddo

  !Mass1
  mass1=0d0
  porLam=0d0  
  do i=1,nelem    !go to each ele
    nq=conn(i,:)  !get node conn
    porCell=(por(nq(1))+por(nq(2))+por(nq(3)))/3d0    
    porCell=((1d0-porCell)**3d0)/porCell    
    do lRow=1,6   !iter thu local row
      gRow=nq(lRow)     !global row
      k=(gRow-1)*ivq(0) !global row pos in 1d global stiff
      nlinkq=linkq(k+1:k+ivq(0))  !link for grow
      do lCol=1,6         !iter thru local col
        gCol=nq(lCol)     !global col
        do j=1,ivq(gRow)  !search for gCol in link of gRow
          if(nlinkq(j).eq.gCol) goto 11
        enddo
        write(9,*)"[Err] node conn missing in Mass1 at",gRow
        stop
        11 mass1(k+j)=mass1(k+j)+(locM1(lRow,lCol)*invJ(i,5))
        porLam(k+j)=porLam(k+j)+(locM1(lRow,lCol)*invJ(i,5) &
          *stoneCnst(5)*porCell)        
      enddo
    enddo
  enddo
  write(9,*) "Mass1 Done"
  ! do i=1,npoint
  !   k=(i-1)*ivq(0)
  !   write(9,*) i,":",mass1(k+1:k+ivq(0))
  ! enddo
  ! write(9,*)

  !Mass2
  mass2=0d0
  do i=1,nelem
    nl=conn(i,1:3)
    do lRow=1,3
      gRow=nl(lRow)
      k=(gRow-1)*ivl(0)
      nlinkl=linkl(k+1:k+ivl(0))
      do lCol=1,3
        gCol=nl(lCol)
        do j=1,ivl(gRow)
          if(nlinkl(j).eq.gCol) goto 12
        enddo
        write(9,*)"[Err] node conn missing in Mass2 at",gRow
        stop
        12 mass2(k+j)=mass2(k+j)+(locM2(lRow,lCol)*invJ(i,5))
      enddo
    enddo
  enddo
  write(9,*) "Mass2 Done"
  ! do i=1,npoinl
  !   k=(i-1)*ivl(0)
  !   write(9,*) i,":",mass2(k+1:k+ivl(0))
  ! enddo
  ! write(9,*)

  !6x3
  !Mass Tou
  touMass=0d0
  do i=1,nelem
    nq=conn(i,:)
    do lRow=1,6
      gRow=nq(lRow)
      k=(gRow-1)*ivl(0)
      nlinkl=linkl(k+1:k+ivl(0))
      do lCol=1,3
        gCol=nq(lCol)
        do j=1,ivl(gRow)
          if(nlinkl(j).eq.gCol) goto 14
        enddo
        write(9,*)"[Err] node conn missing in Mass2 at",gRow
        stop
        14 touMass(k+j)=touMass(k+j) & 
          +(loctouMass(lRow,lCol)*invJ(i,5))
      enddo
    enddo
  enddo
  write(9,*) "Mass2 Done"

end subroutine massMatrices

subroutine matrixSet1(npoinl,npoint,nelem,conn,ivl,ivq,linkl,&
  linkq,invJ,depth,gBs1,gBs2,gBs3,gBs4,gCxFlux,gCyFlux,gAx,gAy,&
  gDMat,por,cnst)
implicit none

  !Bsnq matrices 1,2,3,4
  !c Flux gradient mtrices x,y
  !A matrices x,y
  !D matrix
  
  integer(kind=4),intent(in)::npoinl,npoint,nelem,conn(nelem,6)
  integer(kind=4),intent(in)::ivl(0:npoint),linkl(ivl(0)*npoint)
  integer(kind=4),intent(in)::ivq(0:npoint),linkq(ivq(0)*npoint)
  integer(kind=4)::i,j,k,i2,j2,k2,n(6),gRow,gCol,lRow,lCol
  integer(kind=4)::nlinkl(ivl(0)),nlinkq(ivq(0))

  real(kind=8),intent(in)::invJ(nelem,5),depth(npoint),cnst(10)
  real(kind=8),intent(in)::por(npoinl)
  real(kind=8),intent(out)::gBs1(ivq(0)*npoint)
  real(kind=8),intent(out)::gBs2(ivq(0)*npoint)
  real(kind=8),intent(out)::gBs3(ivq(0)*npoint)
  real(kind=8),intent(out)::gBs4(ivq(0)*npoint)
  real(kind=8),intent(out)::gCxFlux(ivq(0)*npoinl)
  real(kind=8),intent(out)::gCyFlux(ivq(0)*npoinl)
  real(kind=8),intent(out)::gAx(ivl(0)*npoint),gAy(ivl(0)*npoint)
  real(kind=8),intent(out)::gDMat(ivl(0)*npoinl)    
  real(kind=8)::bs1t1(6,6),bs1t2(6,6)
  real(kind=8)::bs2t1(6,6),bs2t2(6,6),bs2t3(6,6)
  real(kind=8)::bs3t1(6,6),bs3t2(6,6),bs3t3(6,6)
  real(kind=8)::bs4t1(6,6),bs4t2(6,6)
  real(kind=8)::cxFlux(3,6),cyFlux(3,6)
  real(kind=8)::ax(6,3),ay(6,3),dMat(3,3),porCell

  gBs1=0d0
  gBs2=0d0
  gBs3=0d0
  gBs4=0d0
  gCxFlux=0d0
  gCyFlux=0d0
  gAx=0d0
  gAy=0d0
  gDMat=0d0

  do i=1,nelem
    n=conn(i,:)    
    call bsnq1Term1(bs1t1,invJ(i,1),invJ(i,2),invJ(i,3),&
      invJ(i,4),depth(n(1)),depth(n(2)),depth(n(3)))

    call bsnq1Term2(bs1t2,invJ(i,1),invJ(i,2),invJ(i,3),&
      invJ(i,4),depth(n(1)),depth(n(2)),depth(n(3)))

    call bsnq2Term1(bs2t1,invJ(i,1),invJ(i,2),invJ(i,3),&
      invJ(i,4),depth(n(1)),depth(n(2)),depth(n(3)))

    call bsnq2Term2(bs2t2,invJ(i,1),invJ(i,2),invJ(i,3),&
      invJ(i,4),depth(n(1)),depth(n(2)),depth(n(3)))

    call bsnq2Term3(bs2t3,invJ(i,1),invJ(i,2),invJ(i,3),&
      invJ(i,4),depth(n(1)),depth(n(2)),depth(n(3)))

    call bsnq3Term1(bs3t1,invJ(i,1),invJ(i,2),invJ(i,3),&
      invJ(i,4),depth(n(1)),depth(n(2)),depth(n(3)))

    bs3t2=bs2t2
    bs3t3=bs2t3

    call bsnq4Term1(bs4t1,invJ(i,1),invJ(i,2),invJ(i,3),&
      invJ(i,4),depth(n(1)),depth(n(2)),depth(n(3)))

    call bsnq4Term2(bs4t2,invJ(i,1),invJ(i,2),invJ(i,3),&
      invJ(i,4),depth(n(1)),depth(n(2)),depth(n(3)))

    call cFluxMat(cxFlux,cyFlux,invJ(i,1),invJ(i,2),&
      invJ(i,3),invJ(i,4))

    call remainMat(ax,ay,dMat,invJ(i,1),invJ(i,2),invJ(i,3),&
      invJ(i,4),depth(n(1)),depth(n(2)),depth(n(3)))    

    
    porCell=(por(n(1))+por(n(2))+por(n(3)))/3.0d0

    bs1t1=invJ(i,5)*((cnst(3)*bs1t1)+(cnst(4)*bs1t2))
    bs2t1=invJ(i,5)*((cnst(3)*bs2t1)+(cnst(5)*bs2t2)-(bs2t3/6d0))
    bs3t1=invJ(i,5)*((cnst(3)*bs3t1)-(bs3t2/6d0)+(cnst(5)*bs3t3))
    bs4t1=invJ(i,5)*((cnst(3)*bs4t1)+(cnst(4)*bs4t2))
    cxFlux=-invJ(i,5)*cxFlux/porCell
    cyFlux=-invJ(i,5)*cyFlux/porCell
    ax=invJ(i,5)*cnst(6)*ax*porCell
    ay=invJ(i,5)*cnst(6)*ay*porCell
    dMat=invJ(i,5)*dMat    

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
        write(9,*)"[Err] node conn missing in Bsnq at",gRow
        stop
        11  gBs1(k+j)=gBs1(k+j)+bs1t1(lRow,lCol)
        gBs2(k+j)=gBs2(k+j)+bs2t1(lRow,lCol)
        gBs3(k+j)=gBs3(k+j)+bs3t1(lRow,lCol)
        gBs4(k+j)=gBs4(k+j)+bs4t1(lRow,lCol)        
      enddo
    enddo

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

    !3x3
    do lRow=1,3
      gRow=n(lRow)
      k=(gRow-1)*ivl(0)
      nlinkl=linkl(k+1:k+ivl(0))
      do lCol=1,3
        gCol=n(lCol)
        do j=1,ivl(gRow)
          if(nlinkl(j).eq.gCol) goto 14
        enddo
        write(9,*)"[Err] node conn missing in dMat at",gRow
        stop
        14 gDMat(k+j)=gDMat(k+j)+dMat(lRow,lCol)        
      enddo
    enddo
  enddo

end subroutine matrixSet1


subroutine bsnq1Term1(bs1t1,b11,b12,b21,b22,h1,h2,h3)
implicit none

  real(kind=8),intent(out)::bs1t1(6,6)
  real(kind=8),intent(in)::b11,b12,b21,b22,h1,h2,h3 

  bs1t1(1,1)=(1d0/180d0)*(b11 + b12)**2d0 *(39d0*h1**2d0 &
    + 15d0*h1*(h2 + h3) + 7d0*(h2**2d0 + h2*h3 + h3**2d0))

  bs1t1(1,2)=(1d0/180d0)*b11*(b11 + b12)*(9d0*h1**2d0 &
    + 9d0*h2**2d0 + 5d0*h2*h3 + h3**2d0 + h1*(h2 + 5d0*h3))

  bs1t1(1,3)=(1d0/180d0)*b12*(b11 + b12)*(9d0*h1**2d0 &
    + h2**2d0 + 5d0*h2*h3 + 9d0*h3**2d0 + h1*(5d0*h2 + h3))

  bs1t1(1,4)=(-(1d0/45d0))*(b11 + b12)*(b11*(12d0*h1**2d0 &
    + 4d0*h1*h2 + 4d0*h2**2d0 + 5d0*h1*h3 + 3d0*h2*h3 &
    + 2d0*h3**2d0) + b12*(-3d0*h1**2d0 + 3d0*h2**2d0 &
    + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)))

  bs1t1(1,5)=(1d0/45d0)*(b11 + b12)*(b12*(-3d0*h1**2d0 &
    + 3d0*h2**2d0 + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)) &
    + b11*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 + 3d0*h3**2d0 &
    - h1*(h2 + 2d0*h3)))

  bs1t1(1,6)=(-(1d0/45d0))*(b11 + b12)*(b12*(12d0*h1**2d0 &
    + 5d0*h1*h2 + 2d0*h2**2d0 + 4d0*h1*h3 + 3d0*h2*h3 &
    + 4d0*h3**2d0) + b11*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 - h1*(h2 + 2d0*h3)))

  bs1t1(2,1)=(1d0/180d0)*b11*(b11 + b12)*(9d0*h1**2d0 &
    + 9d0*h2**2d0 + 5d0*h2*h3 + h3**2d0 + h1*(h2 + 5d0*h3))

  bs1t1(2,2)=(1d0/180d0)*b11**2d0*(7d0*h1**2d0 + 15d0*h1*h2 &
    + 39d0*h2**2d0 + 7d0*h1*h3 + 15d0*h2*h3 + 7d0*h3**2d0)

  bs1t1(2,3)=(-(1d0/180d0))*b11*b12*(h1**2d0 + 9d0*h2**2d0 &
    + h2*h3 + 9d0*h3**2d0 + 5d0*h1*(h2 + h3))

  bs1t1(2,4)=(-(1d0/45d0))*b11*(b11*(4d0*h1**2d0 + 4d0*h1*h2 &
    + 12d0*h2**2d0 + 3d0*h1*h3 + 5d0*h2*h3 + 2d0*h3**2d0) &
    + b12*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs1t1(2,5)=(1d0/45d0)*b11*(b11*(-h1**2d0 + 3d0*h2**2d0 &
    + h1*(h2 - 2d0*h3) + 2d0*h2*h3 - 3d0*h3**2d0) &
    + b12*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs1t1(2,6)=(-(1d0/45d0))*b11*(b12*(h1 - h3)*(2d0*h1 - h2 &
    + 2d0*h3) + b11*(-h1**2d0 + 3d0*h2**2d0 + h1*(h2 - 2d0*h3) &
    + 2d0*h2*h3 - 3d0*h3**2d0))

  bs1t1(3,1)=(1d0/180d0)*b12*(b11 + b12)*(9d0*h1**2d0 + h2**2d0 &
    + 5d0*h2*h3 + 9d0*h3**2d0 + h1*(5d0*h2 + h3))

  bs1t1(3,2)=(-(1d0/180d0))*b11*b12*(h1**2d0 + 9d0*h2**2d0 &
    + h2*h3 + 9d0*h3**2d0 + 5d0*h1*(h2 + h3))

  bs1t1(3,3)=(1d0/180d0)*b12**2d0*(7d0*(h1**2d0 + h1*h2 &
    + h2**2d0) + 15d0*(h1 + h2)*h3 + 39d0*h3**2d0)

  bs1t1(3,4)=(-(1d0/45d0))*b12*(b11*(h1 - h2)*(2d0*h1 + 2d0*h2 &
    - h3) + b12*(-h1**2d0 - 3d0*h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 + h1*(-2d0*h2 + h3)))

  bs1t1(3,5)=(-(1d0/45d0))*b12*(b12*(h1**2d0 + 2d0*h1*h2 &
    + 3d0*h2**2d0 - h1*h3 - 2d0*h2*h3 - 3d0*h3**2d0) &
    - b11*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs1t1(3,6)=(-(1d0/45d0))*b12*(b12*(4d0*h1**2d0 + 3d0*h1*h2 &
    + 2d0*h2**2d0 + 4d0*h1*h3 + 5d0*h2*h3 + 12d0*h3**2d0) &
    + b11*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs1t1(4,1)=(-(1d0/45d0))*(b11 + b12)*(b11*(12d0*h1**2d0 &
    + 4d0*h1*h2 + 4d0*h2**2d0 + 5d0*h1*h3 + 3d0*h2*h3 &
    + 2d0*h3**2d0) + b12*(-3d0*h1**2d0 + 3d0*h2**2d0 &
    + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)))

  bs1t1(4,2)=(-(1d0/45d0))*b11*(b11*(4d0*h1**2d0 + 4d0*h1*h2 &
    + 12d0*h2**2d0 + 3d0*h1*h3 + 5d0*h2*h3 + 2d0*h3**2d0) &
    + b12*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs1t1(4,3)=(-(1d0/45d0))*b12*(b11*(h1 - h2)*(2d0*h1 + 2d0*h2 &
    - h3) + b12*(-h1**2d0 - 3d0*h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 + h1*(-2d0*h2 + h3)))

  bs1t1(4,4)=(4d0/45d0)*(b11*b12*(-h1**2d0 + 2d0*h1*h2 &
    + 9d0*h2**2d0 + 4d0*h2*h3 + h3**2d0) + b11**2d0*(4d0*h1**2d0 &
    + 4d0*h2**2d0 + 2d0*h2*h3 + h3**2d0 + 2d0*h1*(h2 + h3)) &
    + b12**2d0*(h1**2d0 + 6d0*h2**2d0 + 3d0*h2*h3 + h3**2d0 &
    + h1*(3d0*h2 + h3)))

  bs1t1(4,5)=(-(4d0/45d0))*((-b11**2d0)*(h1 - h2)*(h1 + h2 + h3) &
    + b11*b12*(2d0*h1*h2 + 6d0*h2**2d0 + h1*h3 + 4d0*h2*h3 &
    + 2d0*h3**2d0) + b12**2d0*(h1**2d0 + 6d0*h2**2d0 + 3d0*h2*h3 &
    + h3**2d0 + h1*(3d0*h2 + h3)))

  bs1t1(4,6)=(4d0/45d0)*((-b11**2d0)*(h1 - h2)*(h1 + h2 + h3) &
    - b12**2d0*(h1 - h3)*(h1 + h2 + h3) + b11*b12*(4d0*h1**2d0 &
    + 2d0*h2**2d0 + 3d0*h2*h3 + 2d0*h3**2d0 + 2d0*h1*(h2 + h3)))

  bs1t1(5,1)=(1d0/45d0)*(b11 + b12)*(b12*(-3d0*h1**2d0 &
    + 3d0*h2**2d0 + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)) &
    + b11*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 + 3d0*h3**2d0 &
    - h1*(h2 + 2d0*h3)))

  bs1t1(5,2)=(1d0/45d0)*b11*(b11*(-h1**2d0 + 3d0*h2**2d0 &
    + h1*(h2 - 2d0*h3) + 2d0*h2*h3 - 3d0*h3**2d0) &
    + b12*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs1t1(5,3)=(-(1d0/45d0))*b12*(b12*(h1**2d0 + 2d0*h1*h2 &
    + 3d0*h2**2d0 - h1*h3 - 2d0*h2*h3 - 3d0*h3**2d0) &
    - b11*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs1t1(5,4)=(-(4d0/45d0))*((-b11**2d0)*(h1 - h2)*(h1 + h2 + h3) &
    + b11*b12*(2d0*h1*h2 + 6d0*h2**2d0 + h1*h3 + 4d0*h2*h3 &
    + 2d0*h3**2d0) +  b12**2d0*(h1**2d0 + 6d0*h2**2d0 &
    + 3d0*h2*h3 + h3**2d0 + h1*(3d0*h2 + h3)))

  bs1t1(5,5)=(4d0/45d0)*(b11*b12*(h1**2d0 + 3d0*h2**2d0 &
    + 4d0*h2*h3 + 3d0*h3**2d0 + 2d0*h1*(h2 + h3)) &
    + b12**2d0*(h1**2d0 + 6d0*h2**2d0 + 3d0*h2*h3 + h3**2d0 &
    + h1*(3d0*h2 + h3)) + b11**2d0*(h1**2d0 + h2**2d0 &
    + 3d0*h2*h3 + 6d0*h3**2d0 + h1*(h2 + 3d0*h3)))

  bs1t1(5,6)=(-(4d0/45d0))*((-b12**2d0)*(h1 - h3)*(h1 + h2 &
    + h3) + b11*b12*(h1*h2 + 2d0*h2**2d0 + 2d0*h1*h3 &
    + 4d0*h2*h3 + 6d0*h3**2) +  b11**2d0*(h1**2d0 + h2**2d0 &
    + 3d0*h2*h3 + 6d0*h3**2d0 + h1*(h2 + 3d0*h3)))

  bs1t1(6,1)=(-(1d0/45d0))*(b11 + b12)*(b12*(12d0*h1**2d0 &
    + 5d0*h1*h2 + 2d0*h2**2d0 + 4d0*h1*h3 + 3d0*h2*h3 &
    + 4d0*h3**2d0) + b11*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 - h1*(h2 + 2d0*h3)))

  bs1t1(6,2)=(-(1d0/45d0))*b11*(b12*(h1 - h3)*(2d0*h1 - h2 &
    + 2d0*h3) + b11*(-h1**2d0 + 3d0*h2**2d0 + h1*(h2 - 2d0*h3) &
    + 2d0*h2*h3 - 3d0*h3**2d0))

  bs1t1(6,3)=(-(1d0/45d0))*b12*(b12*(4d0*h1**2d0 + 3d0*h1*h2 &
    + 2d0*h2**2d0 + 4d0*h1*h3 + 5d0*h2*h3 + 12d0*h3**2d0) &
    + b11*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs1t1(6,4)=(4d0/45d0)*((-b11**2d0)*(h1 - h2)*(h1 + h2 + h3) &
    - b12**2d0*(h1 - h3)*(h1 + h2 + h3) + b11*b12*(4d0*h1**2d0 &
    + 2d0*h2**2d0 + 3d0*h2*h3 + 2d0*h3**2d0 + 2d0*h1*(h2 + h3)))

  bs1t1(6,5)=(-(4d0/45d0))*((-b12**2d0)*(h1 - h3)*(h1 + h2 + h3) &
    + b11*b12*(h1*h2 + 2d0*h2**2d0 + 2d0*h1*h3 + 4d0*h2*h3 &
    + 6d0*h3**2d0) + b11**2d0*(h1**2d0 + h2**2d0 + 3d0*h2*h3 &
    + 6d0*h3**2d0 + h1*(h2 + 3d0*h3)))

  bs1t1(6,6)=(4d0/45d0)*(b11*b12*(-h1**2d0 + h2**2d0 + 2d0*h1*h3 &
    + 4d0*h2*h3 + 9d0*h3**2) + b12**2d0*(4d0*h1**2d0 + h2**2d0 &
    + 2d0*h2*h3 + 4d0*h3**2d0 + 2d0*h1*(h2 + h3)) &
    + b11**2d0*(h1**2d0 + h2**2d0 + 3d0*h2*h3 + 6d0*h3**2d0 &
    + h1*(h2 + 3d0*h3))) 

end subroutine bsnq1Term1


subroutine bsnq1Term2(bs1t2,b11,b12,b21,b22,h1,h2,h3)
implicit none

  real(kind=8),intent(out)::bs1t2(6,6)
  real(kind=8),intent(in)::b11,b12,b21,b22,h1,h2,h3 

  bs1t2(1,1)=(1d0/120d0)*(b11 + b12)*(b11*(h1 - h2) &
    + b12*(h1 - h3))*(6d0*h1 + h2 + h3)

  bs1t2(1,2)=(1d0/360d0)*b11*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(6d0*h1 + 5d0*h2 + h3)

  bs1t2(1,3)=(1d0/360d0)*b12*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(6d0*h1 + h2 + 5d0*h3)

  bs1t2(1,4)=(-(1d0/90d0))*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b12*(2d0*h2 + h3) + b11*(6d0*h1 + 2d0*h2 + h3))

  bs1t2(1,5)=(1d0/90d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b12*(2d0*h2 + h3) + b11*(h2 + 2d0*h3))

  bs1t2(1,6)=(-(1d0/90d0))*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b11*(h2 + 2d0*h3) + b12*(6d0*h1 + h2 + 2d0*h3))

  bs1t2(2,1)=(1d0/360d0)*(b11 + b12)*(5d0*h1 + 6d0*h2 &
    + h3)*(b11*(-h1 + h2) + b12*(-h1 + h3))

  bs1t2(2,2)=(-(1d0/120d0))*b11*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(h1 + 6d0*h2 + h3)

  bs1t2(2,3)=(1d0/360d0)*b12*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(h1 + 6d0*h2 + 5d0*h3)

  bs1t2(2,4)=(1d0/90d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(6d0*b12*h2 + b11*(2d0*h1 + 6d0*h2 + h3))

  bs1t2(2,5)=(1d0/90d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(-6d0*b12*h2 + b11*(h1 + 2d0*h3))

  bs1t2(2,6)=(-(1d0/90d0))*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b12*(-h1 + h3) + b11*(h1 + 2d0*h3))

  bs1t2(3,1)=(1d0/360d0)*(b11 + b12)*(5d0*h1 + h2 &
    + 6d0*h3)*(b11*(-h1 + h2) + b12*(-h1 + h3))

  bs1t2(3,2)=(1d0/360d0)*b11*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(h1 + 5d0*h2 + 6d0*h3)

  bs1t2(3,3)=(1d0/120d0)*b12*(h1 + h2 + 6d0*h3)*(b11*(-h1 &
    + h2) + b12*(-h1 + h3))

  bs1t2(3,4)=(1d0/90d0)*(b11*(h1 - h2) - b12*(h1 &
    + 2d0*h2))*(b11*(h1 - h2) + b12*(h1 - h3))

  bs1t2(3,5)=(1d0/90d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b12*(h1 + 2d0*h2) - 6d0*b11*h3)

  bs1t2(3,6)=(1d0/90d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(6d0*b11*h3 + b12*(2d0*h1 + h2 + 6d0*h3))

  bs1t2(4,1)=(1d0/90d0)*(b11 + b12)*(b11*(h1 - h2) &
    + b12*(h1 - h3))*(6d0*h1 + 2d0*h2 + h3)

  bs1t2(4,2)=(-(1d0/90d0))*b11*(b11*(h1 - h2) &
    + b12*(h1 - h3))*(2d0*h1 + 6d0*h2 + h3)

  bs1t2(4,3)=(1d0/90d0)*b12*(-2d0*h1 - 2d0*h2 &
    + h3)*(b11*(-h1 + h2) + b12*(-h1 + h3))

  bs1t2(4,4)=(2d0/45d0)*(b11*(-h1 + h2) + b12*(-h1 &
    + h3))*(b11*(h1 - h2) - b12*(2d0*h1 + 3d0*h2 + h3))

  bs1t2(4,5)=(2d0/45d0)*(b11*(-h1 + h2) + b12*(-h1 &
    + h3))*(b11*(h1 + h2 + h3) + b12*(2d0*h1 + 3d0*h2 + h3))

  bs1t2(4,6)=(2d0/45d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*((-b12)*(2d0*h1 + h2) + b11*(h1 + h2 + h3))

  bs1t2(5,1)=(1d0/90d0)*(b11 + b12)*(b11*(h1 - h2) &
    + b12*(h1 - h3))*(h1 - 2d0*(h2 + h3))

  bs1t2(5,2)=(-(1d0/90d0))*b11*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(h1 + 6d0*h2 + 2d0*h3)

  bs1t2(5,3)=(1d0/90d0)*b12*(h1 + 2d0*h2 + 6d0*h3)*(b11*(-h1 &
    + h2) + b12*(-h1 + h3))

  bs1t2(5,4)=(2d0/45d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b11*(2d0*h2 + h3) + b12*(h1 + 3d0*h2 + 2d0*h3))

  bs1t2(5,5)=(2d0/45d0)*(b11*(-h1 + h2) + b12*(-h1 &
    + h3))*(b12*(h1 + 3d0*h2 + 2d0*h3)+b11*(h1 + 2d0*h2 + 3d0*h3))

  bs1t2(5,6)=(2d0/45d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b12*(h2 + 2d0*h3) + b11*(h1 + 2d0*h2 + 3d0*h3))

  bs1t2(6,1)=(1d0/90d0)*(b11 + b12)*(b11*(h1 - h2) &
    + b12*(h1 - h3))*(6d0*h1 + h2 + 2d0*h3)

  bs1t2(6,2)=(1d0/90d0)*b11*(-2d0*h1 + h2 - 2d0*h3)*(b11*(-h1 &
    + h2) + b12*(-h1 + h3))

  bs1t2(6,3)=(-(1d0/90d0))*b12*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(2d0*h1 + h2 + 6d0*h3)

  bs1t2(6,4)=(-(2d0/45d0))*(b11*(-h1 + h2) + b12*(-h1 &
    + h3))*((-b11)*(2d0*h1 + h3) + b12*(h1 + h2 + h3))

  bs1t2(6,5)=(-(2d0/45d0))*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b12*(h1 + h2 + h3) + b11*(2d0*h1 + h2 + 3d0*h3))

  bs1t2(6,6)=(2d0/45d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b12*(-h1 + h3) + b11*(2d0*h1 + h2 + 3d0*h3))

end subroutine bsnq1Term2

subroutine bsnq2Term1(bs2t1,b11,b12,b21,b22,h1,h2,h3)
implicit none

  real(kind=8),intent(out)::bs2t1(6,6)
  real(kind=8),intent(in)::b11,b12,b21,b22,h1,h2,h3 

  bs2t1(1,1)=(1d0/180d0)*(b11 + b12)*(b21 + b22)*(39d0*h1**2d0 &
    + 15d0*h1*(h2 + h3) + 7d0*(h2**2d0 + h2*h3 + h3**2d0))

  bs2t1(1,2)=(1d0/180d0)*(b11 + b12)*b21*(9d0*h1**2d0 &
    + 9d0*h2**2d0 + 5d0*h2*h3 + h3**2d0 + h1*(h2 + 5d0*h3))

  bs2t1(1,3)=(1d0/180d0)*(b11 + b12)*b22*(9d0*h1**2d0 &
    + h2**2d0 + 5d0*h2*h3 + 9d0*h3**2d0 + h1*(5d0*h2 + h3))

  bs2t1(1,4)=(-(1d0/45d0))*(b11 + b12)*(b21*(12d0*h1**2d0 &
    + 4d0*h1*h2 + 4d0*h2**2d0 + 5d0*h1*h3 + 3d0*h2*h3 &
    + 2d0*h3**2d0) + b22*(-3d0*h1**2d0 + 3d0*h2**2d0 &
    + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)))

  bs2t1(1,5)=(1d0/45d0)*(b11 + b12)*(b22*(-3d0*h1**2d0 &
    + 3d0*h2**2d0 + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)) &
    + b21*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 + 3d0*h3**2d0 &
    - h1*(h2 + 2d0*h3)))

  bs2t1(1,6)=(-(1d0/45d0))*(b11 + b12)*(b22*(12d0*h1**2d0 &
    + 5d0*h1*h2 + 2d0*h2**2d0 + 4d0*h1*h3 + 3d0*h2*h3 &
    + 4d0*h3**2d0) + b21*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 - h1*(h2 + 2d0*h3)))

  bs2t1(2,1)=(1d0/180d0)*b11*(b21 + b22)*(9d0*h1**2d0 &
    + 9d0*h2**2d0 + 5d0*h2*h3 + h3**2d0 + h1*(h2 + 5d0*h3))

  bs2t1(2,2)=(1d0/180d0)*b11*b21*(7d0*h1**2d0 + 15d0*h1*h2 &
    + 39d0*h2**2d0 + 7d0*h1*h3 + 15d0*h2*h3 + 7d0*h3**2d0)

  bs2t1(2,3)=(-(1d0/180d0))*b11*b22*(h1**2d0 + 9d0*h2**2d0 &
    + h2*h3 + 9d0*h3**2d0 + 5d0*h1*(h2 + h3))

  bs2t1(2,4)=(-(1d0/45d0))*b11*(b21*(4d0*h1**2d0 + 4d0*h1*h2 &
    + 12d0*h2**2d0 + 3d0*h1*h3 + 5d0*h2*h3 + 2d0*h3**2d0) &
    + b22*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs2t1(2,5)=(1d0/45d0)*b11*(b21*(-h1**2d0 + 3d0*h2**2d0 &
    + h1*(h2 - 2d0*h3) + 2d0*h2*h3 - 3d0*h3**2d0) &
    + b22*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs2t1(2,6)=(-(1d0/45d0))*b11*(b22*(h1 - h3)*(2d0*h1 - h2 &
    + 2d0*h3) + b21*(-h1**2d0 + 3d0*h2**2d0 + h1*(h2 - 2d0*h3) &
    + 2d0*h2*h3 - 3d0*h3**2d0))

  bs2t1(3,1)=(1d0/180d0)*b12*(b21 + b22)*(9d0*h1**2d0 + h2**2d0 &
    + 5d0*h2*h3 + 9d0*h3**2d0 + h1*(5d0*h2 + h3))

  bs2t1(3,2)=(-(1d0/180d0))*b12*b21*(h1**2d0 + 9d0*h2**2d0 &
    + h2*h3 + 9d0*h3**2d0 + 5d0*h1*(h2 + h3))

  bs2t1(3,3)=(1d0/180d0)*b12*b22*(7d0*h1**2d0 + 7d0*h1*h2 &
    + 7d0*h2**2d0 + 15d0*h1*h3 + 15d0*h2*h3 + 39d0*h3**2d0)

  bs2t1(3,4)=(-(1d0/45d0))*b12*(b21*(h1 - h2)*(2d0*h1 + 2d0*h2 &
    - h3) + b22*(-h1**2d0 - 3d0*h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 + h1*(-2d0*h2 + h3)))

  bs2t1(3,5)=(1d0/45d0)*b12*(b22*(-h1**2d0 - 3d0*h2**2d0 &
    + 2d0*h2*h3 + 3d0*h3**2d0 + h1*(-2d0*h2 + h3)) &
    + b21*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs2t1(3,6)=(-(1d0/45d0))*b12*(b22*(4d0*h1**2d0 + 3d0*h1*h2 &
    + 2d0*h2**2d0 + 4d0*h1*h3 + 5d0*h2*h3 + 12d0*h3**2d0) &
    + b21*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs2t1(4,1)=(-(1d0/45d0))*(b21 + b22)*(b11*(12d0*h1**2d0 &
    + 4d0*h1*h2 + 4d0*h2**2d0 + 5d0*h1*h3 + 3d0*h2*h3 &
    + 2d0*h3**2d0) + b12*(-3d0*h1**2d0 + 3d0*h2**2d0 &
    + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)))

  bs2t1(4,2)=(-(1d0/45d0))*b21*(b11*(4d0*h1**2d0 + 4d0*h1*h2 &
    + 12d0*h2**2d0 + 3d0*h1*h3 + 5d0*h2*h3 + 2d0*h3**2d0) &
    + b12*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs2t1(4,3)=(-(1d0/45d0))*b22*(b11*(h1 - h2)*(2d0*h1 + 2d0*h2 &
    - h3) + b12*(-h1**2d0 - 3d0*h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 + h1*(-2d0*h2 + h3)))

  bs2t1(4,4)=(2d0/45d0)*(b12*(2d0*b22*(h1**2d0 + 3d0*h1*h2 &
    + 6d0*h2**2d0 + h1*h3 + 3d0*h2*h3 + h3**2d0) &
    + b21*(-h1**2d0 + 2d0*h1*h2 + 9d0*h2**2d0 + 4d0*h2*h3 &
    + h3**2d0)) + b11*(b22*(-h1**2d0 + 2d0*h1*h2 + 9d0*h2**2d0 &
    + 4d0*h2*h3 + h3**2d0) + 2d0*b21*(4d0*h1**2d0 + 4d0*h2**2d0 &
    + 2d0*h2*h3 + h3**2d0 + 2d0*h1*(h2 + h3))))

  bs2t1(4,5)=(-(2d0/45d0))*(b11*(-2d0*b21*(h1 - h2)*(h1 + h2 &
    + h3) + b22*(-h1**2d0 + 2d0*h1*h2 + 9d0*h2**2d0 + 4d0*h2*h3 &
    + h3**2d0)) + b12*(b21*(h1**2d0 + 3d0*h2**2d0 + 4d0*h2*h3 &
    + 3d0*h3**2d0 + 2d0*h1*(h2 + h3)) + 2d0*b22*(h1**2d0 &
    + 6d0*h2**2d0 + 3d0*h2*h3 + h3**2d0 + h1*(3d0*h2 + h3))))

  bs2t1(4,6)=(2d0/45d0)*(b12*(-2d0*b22*(h1 - h3)*(h1 + h2 + h3) &
    + b21*(h1**2d0 + 3d0*h2**2d0 + 4d0*h2*h3 + 3d0*h3**2d0 &
    + 2d0*h1*(h2 + h3))) + b11*(-2d0*b21*(h1 - h2)*(h1 + h2 &
    + h3) + b22*(7d0*h1**2d0 + 2d0*h1*(h2 + h3) + (h2 + h3)**2d0)))

  bs2t1(5,1)=(1d0/45d0)*(b21 + b22)*(b12*(-3d0*h1**2d0 &
    + 3d0*h2**2d0 + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)) &
    + b11*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 + 3d0*h3**2d0 &
    - h1*(h2 + 2d0*h3)))

  bs2t1(5,2)=(1d0/45d0)*b21*(b11*(-h1**2d0 + 3d0*h2**2d0 &
    + h1*(h2 - 2d0*h3) + 2d0*h2*h3 - 3d0*h3**2d0) + b12*(h1**2d0 &
    + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 + h1*(6d0*h2 + h3)))

  bs2t1(5,3)=(1d0/45d0)*b22*(b12*(-h1**2d0 - 3d0*h2**2d0 &
    + 2d0*h2*h3 + 3d0*h3**2d0 + h1*(-2d0*h2 + h3)) &
    + b11*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs2t1(5,4)=(-(2d0/45d0))*(b12*(2d0*b22*(h1**2d0 + 3d0*h1*h2 &
    + 6d0*h2**2d0 + h1*h3 + 3d0*h2*h3 + h3**2d0) + b21*(-h1**2d0 &
    + 2d0*h1*h2 + 9d0*h2**2d0 + 4d0*h2*h3 + h3**2d0)) &
    + b11*(-2d0*b21*(h1 - h2)*(h1 + h2 + h3) + b22*(h1**2d0 &
    + 3d0*h2**2d0 + 4d0*h2*h3 + 3d0*h3**2d0 + 2d0*h1*(h2 + h3))))

  bs2t1(5,5)=(2d0/45d0)*(b12*(b21*(h1**2d0 + 3d0*h2**2d0 &
    + 4d0*h2*h3 + 3d0*h3**2d0 + 2d0*h1*(h2 + h3)) &
    + 2d0*b22*(h1**2d0 + 6d0*h2**2d0 + 3d0*h2*h3 + h3**2d0 &
    + h1*(3d0*h2 + h3))) + b11*(b22*(h1**2d0 + 3d0*h2**2d0 &
    + 4d0*h2*h3 + 3d0*h3**2d0 + 2d0*h1*(h2 + h3)) &
    + 2d0*b21*(h1**2d0 + h2**2d0 + 3d0*h2*h3 + 6d0*h3**2d0 &
    + h1*(h2 + 3d0*h3))))

  bs2t1(5,6)=(-(2d0/45d0))*(b11*(2d0*b21*(h1**2d0 + h1*h2 &
    + h2**2d0 + 3d0*h1*h3 + 3d0*h2*h3 + 6d0*h3**2d0) &
    + b22*(-h1**2d0 + h2**2d0 + 2d0*h1*h3 + 4d0*h2*h3 &
    + 9d0*h3**2d0)) + b12*(-2d0*b22*(h1 - h3)*(h1 + h2 + h3) &
    + b21*(h1**2d0 + 3d0*h2**2d0 + 4d0*h2*h3 + 3d0*h3**2d0 &
    + 2d0*h1*(h2 + h3))))

  bs2t1(6,1)=(-(1d0/45d0))*(b21 + b22)*(b12*(12d0*h1**2d0 &
    + 5d0*h1*h2 + 2d0*h2**2d0 + 4d0*h1*h3 + 3d0*h2*h3 &
    + 4d0*h3**2d0) + b11*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 - h1*(h2 + 2d0*h3)))

  bs2t1(6,2)=(-(1d0/45d0))*b21*(b12*(h1 - h3)*(2d0*h1 - h2 &
    + 2d0*h3) + b11*(-h1**2d0 + 3d0*h2**2d0 + h1*(h2 - 2d0*h3) &
    + 2d0*h2*h3 - 3d0*h3**2d0))

  bs2t1(6,3)=(-(1d0/45d0))*b22*(b12*(4d0*h1**2d0 + 3d0*h1*h2 &
    + 2d0*h2**2d0 + 4d0*h1*h3 + 5d0*h2*h3 + 12d0*h3**2d0) &
    + b11*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs2t1(6,4)=(2d0/45d0)*(b11*(-2d0*b21*(h1 - h2)*(h1 + h2 + h3) &
    + b22*(h1**2d0 + 3d0*h2**2d0 + 4d0*h2*h3 + 3d0*h3**2d0 &
    + 2d0*h1*(h2 + h3))) + b12*(-2d0*b22*(h1 - h3)*(h1 + h2 &
    + h3) + b21*(7d0*h1**2d0 + 2d0*h1*(h2 + h3) + (h2 + h3)**2d0)))

  bs2t1(6,5)=(-(2d0/45d0))*(b12*(-2d0*b22*(h1 - h3)*(h1 + h2 &
    + h3) + b21*(-h1**2d0 + h2**2d0 + 2d0*h1*h3 + 4d0*h2*h3 &
    + 9d0*h3**2d0)) + b11*(b22*(h1**2d0 + 3d0*h2**2d0 &
    + 4d0*h2*h3 + 3d0*h3**2d0 + 2d0*h1*(h2 + h3)) &
    + 2d0*b21*(h1**2d0 + h2**2d0 + 3d0*h2*h3 + 6d0*h3**2d0 &
    + h1*(h2 + 3d0*h3))))

  bs2t1(6,6)=(2d0/45d0)*(b11*(2d0*b21*(h1**2d0 + h1*h2 + h2**2d0 &
    + 3d0*h1*h3 + 3d0*h2*h3 + 6d0*h3**2d0) + b22*(-h1**2d0 &
    + h2**2d0 + 2d0*h1*h3 + 4d0*h2*h3 + 9d0*h3**2d0)) &
    + b12*(b21*(-h1**2d0 + h2**2d0 + 2d0*h1*h3 + 4d0*h2*h3 &
    + 9d0*h3**2d0) + 2d0*b22*(4d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 &
    + 4d0*h3**2d0 + 2d0*h1*(h2 + h3))))

end subroutine bsnq2Term1

subroutine bsnq2Term2(bs2t2,b11,b12,b21,b22,h1,h2,h3)
implicit none

  real(kind=8),intent(out)::bs2t2(6,6)
  real(kind=8),intent(in)::b11,b12,b21,b22,h1,h2,h3 

  bs2t2(1,1)=(1d0/120d0)*(b21 + b22)*(b11*(h1 - h2) &
    + b12*(h1 - h3))*(6d0*h1 + h2 + h3)

  bs2t2(1,2)=(1d0/360d0)*b21*(b11*(h1 - h2) &
    + b12*(h1 - h3))*(6d0*h1 + 5d0*h2 + h3)

  bs2t2(1,3)=(1d0/360d0)*b22*(b11*(h1 - h2) &
    + b12*(h1 - h3))*(6d0*h1 + h2 + 5d0*h3)

  bs2t2(1,4)=(-(1d0/90d0))*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b22*(2d0*h2 + h3) + b21*(6d0*h1 + 2d0*h2 + h3))

  bs2t2(1,5)=(1d0/90d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b22*(2d0*h2 + h3) + b21*(h2 + 2d0*h3))

  bs2t2(1,6)=(-(1d0/90d0))*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b21*(h2 + 2d0*h3) + b22*(6d0*h1 + h2 + 2d0*h3))

  bs2t2(2,1)=(1d0/360d0)*(b21 + b22)*(5d0*h1 + 6d0*h2 &
    + h3)*(b11*(-h1 + h2) + b12*(-h1 + h3))

  bs2t2(2,2)=(-(1d0/120d0))*b21*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(h1 + 6d0*h2 + h3)

  bs2t2(2,3)=(1d0/360d0)*b22*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(h1 + 6d0*h2 + 5d0*h3)

  bs2t2(2,4)=(1d0/90d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(6d0*b22*h2 + b21*(2d0*h1 + 6d0*h2 + h3))

  bs2t2(2,5)=(1d0/90d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(-6d0*b22*h2 + b21*(h1 + 2d0*h3))

  bs2t2(2,6)=(-(1d0/90d0))*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b22*(-h1 + h3) + b21*(h1 + 2d0*h3))

  bs2t2(3,1)=(1d0/360d0)*(b21 + b22)*(5d0*h1 + h2 &
    + 6d0*h3)*(b11*(-h1 + h2) + b12*(-h1 + h3))

  bs2t2(3,2)=(1d0/360d0)*b21*(b11*(h1 - h2) &
    + b12*(h1 - h3))*(h1 + 5d0*h2 + 6d0*h3)

  bs2t2(3,3)=(1d0/120d0)*b22*(h1 + h2 + 6d0*h3)*(b11*(-h1 &
    + h2) + b12*(-h1 + h3))

  bs2t2(3,4)=(1d0/90d0)*(b21*(h1 - h2) - b22*(h1 &
    + 2d0*h2))*(b11*(h1 - h2) + b12*(h1 - h3))

  bs2t2(3,5)=(1d0/90d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b22*(h1 + 2d0*h2) - 6d0*b21*h3)

  bs2t2(3,6)=(1d0/90d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(6d0*b21*h3 + b22*(2d0*h1 + h2 + 6d0*h3))

  bs2t2(4,1)=(1d0/90d0)*(b21 + b22)*(b11*(h1 - h2) &
    + b12*(h1 - h3))*(6d0*h1 + 2d0*h2 + h3)

  bs2t2(4,2)=(-(1d0/90d0))*b21*(b11*(h1 - h2) &
    + b12*(h1 - h3))*(2d0*h1 + 6d0*h2 + h3)

  bs2t2(4,3)=(1d0/90d0)*b22*(-2d0*h1 - 2d0*h2 + h3)*(b11*(-h1 &
    + h2) + b12*(-h1 + h3))

  bs2t2(4,4)=(2d0/45d0)*(b11*(-h1 + h2) + b12*(-h1 &
    + h3))*(b21*(h1 - h2) - b22*(2d0*h1 + 3d0*h2 + h3))

  bs2t2(4,5)=(2d0/45d0)*(b11*(-h1 + h2) + b12*(-h1 &
    + h3))*(b21*(h1 + h2 + h3) + b22*(2d0*h1 + 3d0*h2 + h3))

  bs2t2(4,6)=(2d0/45d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*((-b22)*(2d0*h1 + h2) + b21*(h1 + h2 + h3))

  bs2t2(5,1)=(1d0/90d0)*(b21 + b22)*(b11*(h1 - h2) &
    + b12*(h1 - h3))*(h1 - 2d0*(h2 + h3))

  bs2t2(5,2)=(-(1d0/90d0))*b21*(b11*(h1 - h2) &
    + b12*(h1 - h3))*(h1 + 6d0*h2 + 2d0*h3)

  bs2t2(5,3)=(1d0/90d0)*b22*(h1 + 2d0*h2 + 6d0*h3)*(b11*(-h1 &
    + h2) + b12*(-h1 + h3))

  bs2t2(5,4)=(2d0/45d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b21*(2d0*h2 + h3) + b22*(h1 + 3d0*h2 + 2d0*h3))

  bs2t2(5,5)=(2d0/45d0)*(b11*(-h1 + h2) + b12*(-h1 &
    + h3))*(b22*(h1 + 3d0*h2 + 2d0*h3) + b21*(h1 + 2d0*h2 &
    + 3d0*h3))

  bs2t2(5,6)=(2d0/45d0)*(b11*(h1 - h2) + b12*(h1 - h3))*(b22*(h2 &
    + 2d0*h3) + b21*(h1 + 2d0*h2 + 3d0*h3))

  bs2t2(6,1)=(1d0/90d0)*(b21 + b22)*(b11*(h1 - h2) &
    + b12*(h1 - h3))*(6d0*h1 + h2 + 2d0*h3)

  bs2t2(6,2)=(1d0/90d0)*b21*(-2d0*h1 + h2 - 2d0*h3)*(b11*(-h1 &
    + h2) + b12*(-h1 + h3))

  bs2t2(6,3)=(-(1d0/90d0))*b22*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(2d0*h1 + h2 + 6d0*h3)

  bs2t2(6,4)=(-(2d0/45d0))*(b11*(-h1 + h2) + b12*(-h1 &
    + h3))*((-b21)*(2d0*h1 + h3) + b22*(h1 + h2 + h3))

  bs2t2(6,5)=(-(2d0/45d0))*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b22*(h1 + h2 + h3) + b21*(2d0*h1 + h2 + 3d0*h3))

  bs2t2(6,6)=(2d0/45d0)*(b11*(h1 - h2) + b12*(h1 &
    - h3))*(b22*(-h1 + h3) + b21*(2d0*h1 + h2 + 3d0*h3))

end subroutine bsnq2Term2

subroutine bsnq2Term3(bs2t3,b11,b12,b21,b22,h1,h2,h3)
implicit none

  real(kind=8),intent(out)::bs2t3(6,6)
  real(kind=8),intent(in)::b11,b12,b21,b22,h1,h2,h3 

  bs2t3(1,1)=(1d0/120d0)*(b11 + b12)*(b21*(h1 - h2) &
    + b22*(h1 - h3))*(6d0*h1 + h2 + h3)

  bs2t3(1,2)=(1d0/360d0)*b11*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(6d0*h1 + 5d0*h2 + h3)

  bs2t3(1,3)=(1d0/360d0)*b12*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(6d0*h1 + h2 + 5d0*h3)

  bs2t3(1,4)=(-(1d0/90d0))*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b12*(2d0*h2 + h3) + b11*(6d0*h1 + 2d0*h2 + h3))

  bs2t3(1,5)=(1d0/90d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b12*(2d0*h2 + h3) + b11*(h2 + 2d0*h3))

  bs2t3(1,6)=(-(1d0/90d0))*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b11*(h2 + 2d0*h3) + b12*(6d0*h1 + h2 + 2d0*h3))

  bs2t3(2,1)=(1d0/360d0)*(b11 + b12)*(5d0*h1 + 6d0*h2 &
    + h3)*(b21*(-h1 + h2) + b22*(-h1 + h3))

  bs2t3(2,2)=(-(1d0/120d0))*b11*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(h1 + 6d0*h2 + h3)

  bs2t3(2,3)=(1d0/360d0)*b12*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(h1 + 6d0*h2 + 5d0*h3)

  bs2t3(2,4)=(1d0/90d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(6d0*b12*h2 + b11*(2d0*h1 + 6d0*h2 + h3))

  bs2t3(2,5)=(1d0/90d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(-6d0*b12*h2 + b11*(h1 + 2d0*h3))

  bs2t3(2,6)=(-(1d0/90d0))*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b12*(-h1 + h3) + b11*(h1 + 2d0*h3))

  bs2t3(3,1)=(1d0/360d0)*(b11 + b12)*(5d0*h1 + h2 &
    + 6d0*h3)*(b21*(-h1 + h2) + b22*(-h1 + h3))

  bs2t3(3,2)=(1d0/360d0)*b11*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(h1 + 5d0*h2 + 6d0*h3)

  bs2t3(3,3)=(1d0/120d0)*b12*(h1 + h2 + 6d0*h3)*(b21*(-h1 &
    + h2) + b22*(-h1 + h3))

  bs2t3(3,4)=(1d0/90d0)*(b11*(h1 - h2) - b12*(h1 &
    + 2d0*h2))*(b21*(h1 - h2) + b22*(h1 - h3))

  bs2t3(3,5)=(1d0/90d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b12*(h1 + 2d0*h2) - 6d0*b11*h3)

  bs2t3(3,6)=(1d0/90d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(6d0*b11*h3 + b12*(2d0*h1 + h2 + 6d0*h3))

  bs2t3(4,1)=(1d0/90d0)*(b11 + b12)*(b21*(h1 - h2) &
    + b22*(h1 - h3))*(6d0*h1 + 2d0*h2 + h3)

  bs2t3(4,2)=(-(1d0/90d0))*b11*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(2d0*h1 + 6d0*h2 + h3)

  bs2t3(4,3)=(1d0/90d0)*b12*(-2d0*h1 - 2d0*h2 + h3)*(b21*(-h1 &
    + h2) + b22*(-h1 + h3))

  bs2t3(4,4)=(2d0/45d0)*(b21*(-h1 + h2) + b22*(-h1 &
    + h3))*(b11*(h1 - h2) - b12*(2d0*h1 + 3d0*h2 + h3))

  bs2t3(4,5)=(2d0/45d0)*(b21*(-h1 + h2) + b22*(-h1 &
    + h3))*(b11*(h1 + h2 + h3) + b12*(2d0*h1 + 3d0*h2 + h3))

  bs2t3(4,6)=(2d0/45d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*((-b12)*(2d0*h1 + h2) + b11*(h1 + h2 + h3))

  bs2t3(5,1)=(1d0/90d0)*(b11 + b12)*(b21*(h1 - h2) &
    + b22*(h1 - h3))*(h1 - 2d0*(h2 + h3))

  bs2t3(5,2)=(-(1d0/90d0))*b11*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(h1 + 6d0*h2 + 2d0*h3)

  bs2t3(5,3)=(1d0/90d0)*b12*(h1 + 2d0*h2 + 6d0*h3)*(b21*(-h1 &
    + h2) + b22*(-h1 + h3))

  bs2t3(5,4)=(2d0/45d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b11*(2d0*h2 + h3) + b12*(h1 + 3d0*h2 + 2d0*h3))

  bs2t3(5,5)=(2d0/45d0)*(b21*(-h1 + h2) + b22*(-h1 &
    + h3))*(b12*(h1 + 3d0*h2 + 2d0*h3) + b11*(h1 + 2d0*h2 &
    + 3d0*h3))

  bs2t3(5,6)=(2d0/45d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b12*(h2 + 2d0*h3) + b11*(h1 + 2d0*h2 + 3d0*h3))

  bs2t3(6,1)=(1d0/90d0)*(b11 + b12)*(b21*(h1 - h2) &
    + b22*(h1 - h3))*(6d0*h1 + h2 + 2d0*h3)

  bs2t3(6,2)=(1d0/90d0)*b11*(-2d0*h1 + h2 - 2d0*h3)*(b21*(-h1 &
    + h2) + b22*(-h1 + h3))

  bs2t3(6,3)=(-(1d0/90d0))*b12*(b21*(h1 - h2) &
    + b22*(h1 - h3))*(2d0*h1 + h2 + 6d0*h3)

  bs2t3(6,4)=(-(2d0/45d0))*(b21*(-h1 + h2) + b22*(-h1 &
    + h3))*((-b11)*(2d0*h1 + h3) + b12*(h1 + h2 + h3))

  bs2t3(6,5)=(-(2d0/45d0))*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b12*(h1 + h2 + h3) + b11*(2d0*h1 + h2 + 3d0*h3))

  bs2t3(6,6)=(2d0/45d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b12*(-h1 + h3) + b11*(2d0*h1 + h2 + 3d0*h3))

end subroutine bsnq2Term3

subroutine bsnq3Term1(bs3t1,b11,b12,b21,b22,h1,h2,h3)
implicit none

  real(kind=8),intent(out)::bs3t1(6,6)
  real(kind=8),intent(in)::b11,b12,b21,b22,h1,h2,h3

  bs3t1(1,1)=(1d0/180d0)*(b11 + b12)*(b21 + b22)*(39d0*h1**2d0 &
    + 15d0*h1*(h2 + h3) + 7d0*(h2**2d0 + h2*h3 + h3**2d0))

  bs3t1(1,2)=(1d0/180d0)*b11*(b21 + b22)*(9d0*h1**2d0 &
    + 9d0*h2**2d0 + 5d0*h2*h3 + h3**2d0 + h1*(h2 + 5d0*h3))

  bs3t1(1,3)=(1d0/180d0)*b12*(b21 + b22)*(9d0*h1**2d0 + h2**2d0 &
    + 5d0*h2*h3 + 9d0*h3**2d0 + h1*(5d0*h2 + h3))

  bs3t1(1,4)=(-(1d0/45d0))*(b21 + b22)*(b11*(12d0*h1**2d0 &
    + 4d0*h1*h2 + 4d0*h2**2d0 + 5d0*h1*h3 + 3d0*h2*h3 &
    + 2d0*h3**2d0) + b12*(-3d0*h1**2d0 + 3d0*h2**2d0 &
    + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)))

  bs3t1(1,5)=(1d0/45d0)*(b21 + b22)*(b12*(-3d0*h1**2d0 &
    + 3d0*h2**2d0 + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 &
    + h3)) + b11*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 - h1*(h2 + 2d0*h3)))

  bs3t1(1,6)=(-(1d0/45d0))*(b21 + b22)*(b12*(12d0*h1**2d0 &
    + 5d0*h1*h2 + 2d0*h2**2d0 + 4d0*h1*h3 + 3d0*h2*h3 &
    + 4d0*h3**2d0) + b11*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 - h1*(h2 + 2d0*h3)))

  bs3t1(2,1)=(1d0/180d0)*(b11 + b12)*b21*(9d0*h1**2d0 &
    + 9d0*h2**2d0 + 5d0*h2*h3 + h3**2d0 + h1*(h2 + 5d0*h3))

  bs3t1(2,2)=(1d0/180d0)*b11*b21*(7d0*h1**2d0 + 15d0*h1*h2 &
    + 39d0*h2**2d0 + 7d0*h1*h3 + 15d0*h2*h3 + 7d0*h3**2d0)

  bs3t1(2,3)=(-(1d0/180d0))*b12*b21*(h1**2d0 + 9d0*h2**2d0 &
    + h2*h3 + 9d0*h3**2d0 + 5d0*h1*(h2 + h3))

  bs3t1(2,4)=(-(1d0/45d0))*b21*(b11*(4d0*h1**2d0 + 4d0*h1*h2 &
    + 12d0*h2**2d0 + 3d0*h1*h3 + 5d0*h2*h3 + 2d0*h3**2d0) &
    + b12*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs3t1(2,5)=(1d0/45d0)*b21*(b11*(-h1**2d0 + 3d0*h2**2d0 &
    + h1*(h2 - 2d0*h3) + 2d0*h2*h3 - 3d0*h3**2d0) &
    + b12*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs3t1(2,6)=(-(1d0/45d0))*b21*(b12*(h1 - h3)*(2d0*h1 &
    - h2 + 2d0*h3) + b11*(-h1**2d0 + 3d0*h2**2d0 + h1*(h2 &
    - 2d0*h3) + 2d0*h2*h3 - 3d0*h3**2d0))

  bs3t1(3,1)=(1d0/180d0)*(b11 + b12)*b22*(9d0*h1**2d0 &
    + h2**2d0 + 5d0*h2*h3 + 9d0*h3**2d0 + h1*(5d0*h2 + h3))

  bs3t1(3,2)=(-(1d0/180d0))*b11*b22*(h1**2d0 + 9d0*h2**2d0 &
    + h2*h3 + 9d0*h3**2d0 + 5d0*h1*(h2 + h3))

  bs3t1(3,3)=(1d0/180d0)*b12*b22*(7d0*h1**2d0 + 7d0*h1*h2 &
    + 7d0*h2**2d0 + 15d0*h1*h3 + 15d0*h2*h3 + 39d0*h3**2d0)

  bs3t1(3,4)=(-(1d0/45d0))*b22*(b11*(h1 - h2)*(2d0*h1 + 2d0*h2 &
    - h3) + b12*(-h1**2d0 - 3d0*h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 + h1*(-2d0*h2 + h3)))

  bs3t1(3,5)=(1d0/45d0)*b22*(b12*(-h1**2d0 - 3d0*h2**2d0 &
    + 2d0*h2*h3 + 3d0*h3**2d0 + h1*(-2d0*h2 + h3)) &
    + b11*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs3t1(3,6)=(-(1d0/45d0))*b22*(b12*(4d0*h1**2d0 + 3d0*h1*h2 &
    + 2d0*h2**2d0 + 4d0*h1*h3 + 5d0*h2*h3 + 12d0*h3**2d0) &
    + b11*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs3t1(4,1)=(-(1d0/45d0))*(b11 + b12)*(b21*(12d0*h1**2d0 &
    + 4d0*h1*h2 + 4d0*h2**2d0 + 5d0*h1*h3 + 3d0*h2*h3 &
    + 2d0*h3**2d0) + b22*(-3d0*h1**2d0 + 3d0*h2**2d0 &
    + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)))

  bs3t1(4,2)=(-(1d0/45d0))*b11*(b21*(4d0*h1**2d0 + 4d0*h1*h2 &
    + 12d0*h2**2d0 + 3d0*h1*h3 + 5d0*h2*h3 + 2d0*h3**2d0) &
    + b22*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs3t1(4,3)=(-(1d0/45d0))*b12*(b21*(h1 - h2)*(2d0*h1 + 2d0*h2 &
    - h3) + b22*(-h1**2d0 - 3d0*h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 + h1*(-2d0*h2 + h3)))

  bs3t1(4,4)=(2d0/45d0)*(b12*(2d0*b22*(h1**2d0 + 3d0*h1*h2 &
    + 6d0*h2**2d0 + h1*h3 + 3d0*h2*h3 + h3**2d0) &
    + b21*(-h1**2d0 + 2d0*h1*h2 + 9d0*h2**2d0 + 4d0*h2*h3 &
    + h3**2d0)) + b11*(b22*(-h1**2d0 + 2d0*h1*h2 + 9d0*h2**2d0 &
    + 4d0*h2*h3 + h3**2d0) + 2d0*b21*(4d0*h1**2d0 + 4d0*h2**2d0 &
    + 2d0*h2*h3 + h3**2d0 + 2d0*h1*(h2 + h3))))

  bs3t1(4,5)=(-(2d0/45d0))*(b12*(2d0*b22*(h1**2d0 + 3d0*h1*h2 &
    + 6d0*h2**2d0 + h1*h3 + 3d0*h2*h3 + h3**2d0) &
    + b21*(-h1**2d0 + 2d0*h1*h2 + 9d0*h2**2d0 + 4d0*h2*h3 &
    + h3**2d0)) + b11*(-2d0*b21*(h1 - h2)*(h1 + h2 + h3) &
    + b22*(h1**2d0 + 3d0*h2**2d0 + 4d0*h2*h3 + 3d0*h3**2d0 &
    + 2d0*h1*(h2 + h3))))

  bs3t1(4,6)=(2d0/45d0)*(b11*(-2d0*b21*(h1 - h2)*(h1 + h2 &
    + h3) + b22*(h1**2d0 + 3d0*h2**2d0 + 4d0*h2*h3 &
    + 3d0*h3**2d0 + 2d0*h1*(h2 + h3))) + b12*(-2d0*b22*(h1 &
    - h3)*(h1 + h2 + h3) + b21*(7d0*h1**2d0 + 2d0*h1*(h2 + h3) &
    + (h2 + h3)**2d0)))

  bs3t1(5,1)=(1d0/45d0)*(b11 + b12)*(b22*(-3d0*h1**2d0 &
    + 3d0*h2**2d0 + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)) &
    + b21*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 + 3d0*h3**2d0 &
    - h1*(h2 + 2d0*h3)))

  bs3t1(5,2)=(1d0/45d0)*b11*(b21*(-h1**2d0 + 3d0*h2**2d0 &
    + h1*(h2 - 2d0*h3) + 2d0*h2*h3 - 3d0*h3**2d0) &
    + b22*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs3t1(5,3)=(1d0/45d0)*b12*(b22*(-h1**2d0 - 3d0*h2**2d0 &
    + 2d0*h2*h3 + 3d0*h3**2d0 + h1*(-2d0*h2 + h3)) &
    + b21*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs3t1(5,4)=(-(2d0/45d0))*(b11*(-2d0*b21*(h1 - h2)*(h1 + h2 &
    + h3) + b22*(-h1**2d0 + 2d0*h1*h2 + 9d0*h2**2d0 + 4d0*h2*h3 &
    + h3**2d0)) + b12*(b21*(h1**2d0 + 3d0*h2**2d0 + 4d0*h2*h3 &
    + 3d0*h3**2d0 + 2d0*h1*(h2 + h3)) + 2d0*b22*(h1**2d0 &
    + 6d0*h2**2d0 + 3d0*h2*h3 + h3**2d0 + h1*(3d0*h2 + h3))))

  bs3t1(5,5)=(2d0/45d0)*(b12*(b21*(h1**2d0 + 3d0*h2**2d0 &
    + 4d0*h2*h3 + 3d0*h3**2d0 + 2d0*h1*(h2 + h3)) &
    + 2d0*b22*(h1**2d0 + 6d0*h2**2d0 + 3d0*h2*h3 + h3**2d0 &
    + h1*(3d0*h2 + h3))) + b11*(b22*(h1**2d0 + 3d0*h2**2d0 &
    + 4d0*h2*h3 + 3d0*h3**2d0 + 2d0*h1*(h2 + h3)) &
    + 2d0*b21*(h1**2d0 + h2**2d0 + 3d0*h2*h3 + 6d0*h3**2d0 &
    + h1*(h2 + 3d0*h3))))

  bs3t1(5,6)=(-(2d0/45d0))*(b12*(-2d0*b22*(h1 - h3)*(h1 + h2 &
    + h3) + b21*(-h1**2d0 + h2**2d0 + 2d0*h1*h3 + 4d0*h2*h3 &
    + 9d0*h3**2d0)) + b11*(b22*(h1**2d0 + 3d0*h2**2d0 &
    + 4d0*h2*h3 + 3d0*h3**2d0 + 2d0*h1*(h2 + h3)) &
    + 2d0*b21*(h1**2d0 + h2**2d0 + 3d0*h2*h3 + 6d0*h3**2d0 &
    + h1*(h2 + 3d0*h3))))

  bs3t1(6,1)=(-(1d0/45d0))*(b11 + b12)*(b22*(12d0*h1**2d0 &
    + 5d0*h1*h2 + 2d0*h2**2d0 + 4d0*h1*h3 + 3d0*h2*h3 &
    + 4d0*h3**2d0) + b21*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 - h1*(h2 + 2d0*h3)))

  bs3t1(6,2)=(-(1d0/45d0))*b11*(b22*(h1 - h3)*(2d0*h1 - h2 &
    + 2d0*h3) + b21*(-h1**2d0 + 3d0*h2**2d0 + h1*(h2 - 2d0*h3) &
    + 2d0*h2*h3 - 3d0*h3**2d0))

  bs3t1(6,3)=(-(1d0/45d0))*b12*(b22*(4d0*h1**2d0 + 3d0*h1*h2 &
    + 2d0*h2**2d0 + 4d0*h1*h3 + 5d0*h2*h3 + 12d0*h3**2d0) &
    + b21*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs3t1(6,4)=(2d0/45d0)*(b12*(-2d0*b22*(h1 - h3)*(h1 + h2 &
    + h3) + b21*(h1**2d0 + 3d0*h2**2d0 + 4d0*h2*h3 + 3d0*h3**2d0 &
    + 2d0*h1*(h2 + h3))) + b11*(-2d0*b21*(h1 - h2)*(h1 + h2 &
    + h3) + b22*(7d0*h1**2d0 + 2d0*h1*(h2 + h3) &
    + (h2 + h3)**2d0)))

  bs3t1(6,5)=(-(2d0/45d0))*(b11*(2d0*b21*(h1**2d0 + h1*h2 &
    + h2**2d0 + 3d0*h1*h3 + 3d0*h2*h3 + 6d0*h3**2d0) &
    + b22*(-h1**2d0 + h2**2d0 + 2d0*h1*h3 + 4d0*h2*h3 &
    + 9d0*h3**2d0)) + b12*(-2d0*b22*(h1 - h3)*(h1 + h2 + h3) &
    + b21*(h1**2d0 + 3d0*h2**2d0 + 4d0*h2*h3 + 3d0*h3**2d0 &
    + 2d0*h1*(h2 + h3))))

  bs3t1(6,6)=(2d0/45d0)*(b11*(2d0*b21*(h1**2d0 + h1*h2 + h2**2d0 &
    + 3d0*h1*h3 + 3d0*h2*h3 + 6d0*h3**2d0) + b22*(-h1**2d0 &
    + h2**2d0 + 2d0*h1*h3 + 4d0*h2*h3 + 9d0*h3**2d0)) &
    + b12*(b21*(-h1**2d0 + h2**2d0 + 2d0*h1*h3 + 4d0*h2*h3 &
    + 9d0*h3**2d0) + 2d0*b22*(4d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 &
    + 4d0*h3**2d0 + 2d0*h1*(h2 + h3))))

end subroutine bsnq3Term1

subroutine bsnq4Term1(bs4t1,b11,b12,b21,b22,h1,h2,h3)
implicit none

  real(kind=8),intent(out)::bs4t1(6,6)
  real(kind=8),intent(in)::b11,b12,b21,b22,h1,h2,h3 

  bs4t1(1,1)=(1d0/180d0)*(b21 + b22)**2d0*(39d0*h1**2d0 &
    + 15d0*h1*(h2 + h3) + 7d0*(h2**2d0 + h2*h3 + h3**2d0))

  bs4t1(1,2)=(1d0/180d0)*b21*(b21 + b22)*(9d0*h1**2d0 &
    + 9d0*h2**2d0 + 5d0*h2*h3 + h3**2d0 + h1*(h2 + 5d0*h3))

  bs4t1(1,3)=(1d0/180d0)*b22*(b21 + b22)*(9d0*h1**2d0 + h2**2d0 &
    + 5d0*h2*h3 + 9d0*h3**2d0 + h1*(5d0*h2 + h3))

  bs4t1(1,4)=(-(1d0/45d0))*(b21 + b22)*(b21*(12d0*h1**2d0 &
    + 4d0*h1*h2 + 4d0*h2**2d0 + 5d0*h1*h3 + 3d0*h2*h3 &
    + 2d0*h3**2d0) + b22*(-3d0*h1**2d0 + 3d0*h2**2d0 &
    + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)))

  bs4t1(1,5)=(1d0/45d0)*(b21 + b22)*(b22*(-3d0*h1**2d0 &
    + 3d0*h2**2d0 + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)) &
    + b21*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 + 3d0*h3**2d0 &
    - h1*(h2 + 2d0*h3)))

  bs4t1(1,6)=(-(1d0/45d0))*(b21 + b22)*(b22*(12d0*h1**2d0 &
    + 5d0*h1*h2 + 2d0*h2**2d0 + 4d0*h1*h3 + 3d0*h2*h3 &
    + 4d0*h3**2d0) + b21*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 - h1*(h2 + 2d0*h3)))

  bs4t1(2,1)=(1d0/180d0)*b21*(b21 + b22)*(9d0*h1**2d0 &
    + 9d0*h2**2d0 + 5d0*h2*h3 + h3**2d0 + h1*(h2 + 5d0*h3))

  bs4t1(2,2)=(1d0/180d0)*b21**2d0*(7d0*h1**2d0 + 15d0*h1*h2 &
    + 39d0*h2**2d0 + 7d0*h1*h3 + 15d0*h2*h3 + 7d0*h3**2d0)

  bs4t1(2,3)=(-(1d0/180d0))*b21*b22*(h1**2d0 + 9d0*h2**2d0 &
    + h2*h3 + 9d0*h3**2d0 + 5d0*h1*(h2 + h3))

  bs4t1(2,4)=(-(1d0/45d0))*b21*(b21*(4d0*h1**2d0 + 4d0*h1*h2 &
    + 12d0*h2**2d0 + 3d0*h1*h3 + 5d0*h2*h3 + 2d0*h3**2d0) &
    + b22*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs4t1(2,5)=(1d0/45d0)*b21*(b21*(-h1**2d0 + 3d0*h2**2d0 &
    + h1*(h2 - 2d0*h3) + 2d0*h2*h3 - 3d0*h3**2d0) &
    + b22*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs4t1(2,6)=(-(1d0/45d0))*b21*(b22*(h1 - h3)*(2d0*h1 - h2 &
    + 2d0*h3) + b21*(-h1**2d0 + 3d0*h2**2d0 + h1*(h2 - 2d0*h3) &
    + 2d0*h2*h3 - 3d0*h3**2d0))

  bs4t1(3,1)=(1d0/180d0)*b22*(b21 + b22)*(9d0*h1**2d0 + h2**2d0 &
    + 5d0*h2*h3 + 9d0*h3**2d0 + h1*(5d0*h2 + h3))

  bs4t1(3,2)=(-(1d0/180d0))*b21*b22*(h1**2d0 + 9d0*h2**2d0 &
    + h2*h3 + 9d0*h3**2d0 + 5d0*h1*(h2 + h3))

  bs4t1(3,3)=(1d0/180d0)*b22**2d0*(7d0*(h1**2d0 + h1*h2 &
    + h2**2d0) + 15d0*(h1 + h2)*h3 + 39d0*h3**2d0)

  bs4t1(3,4)=(-(1d0/45d0))*b22*(b21*(h1 - h2)*(2d0*h1 + 2d0*h2 &
    - h3) + b22*(-h1**2d0 - 3d0*h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 + h1*(-2d0*h2 + h3)))

  bs4t1(3,5)=(-(1d0/45d0))*b22*(b22*(h1**2d0 + 2d0*h1*h2 &
    + 3d0*h2**2d0 - h1*h3 - 2d0*h2*h3 - 3d0*h3**2d0) &
    - b21*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs4t1(3,6)=(-(1d0/45d0))*b22*(b22*(4d0*h1**2d0 + 3d0*h1*h2 &
    + 2d0*h2**2d0 + 4d0*h1*h3 + 5d0*h2*h3 + 12d0*h3**2d0) &
    + b21*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs4t1(4,1)=(-(1d0/45d0))*(b21 + b22)*(b21*(12d0*h1**2d0 &
    + 4d0*h1*h2 + 4d0*h2**2d0 + 5d0*h1*h3 + 3d0*h2*h3 &
    + 2d0*h3**2d0) + b22*(-3d0*h1**2d0 + 3d0*h2**2d0 &
    + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)))

  bs4t1(4,2)=(-(1d0/45d0))*b21*(b21*(4d0*h1**2d0 + 4d0*h1*h2 &
    + 12d0*h2**2d0 + 3d0*h1*h3 + 5d0*h2*h3 + 2d0*h3**2d0) &
    + b22*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs4t1(4,3)=(-(1d0/45d0))*b22*(b21*(h1 - h2)*(2d0*h1 + 2d0*h2 &
    - h3) + b22*(-h1**2d0 - 3d0*h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 + h1*(-2d0*h2 + h3)))

  bs4t1(4,4)=(4d0/45d0)*(b21*b22*(-h1**2d0 + 2d0*h1*h2 &
    + 9d0*h2**2d0 + 4d0*h2*h3 + h3**2d0) &
    + b21**2d0*(4d0*h1**2d0 + 4d0*h2**2d0 + 2d0*h2*h3 &
    + h3**2d0 + 2d0*h1*(h2 + h3)) + b22**2d0*(h1**2d0 &
    + 6d0*h2**2d0 + 3d0*h2*h3 + h3**2d0 + h1*(3d0*h2 + h3)))

  bs4t1(4,5)=(-(4d0/45d0))*((-b21**2d0)*(h1 - h2)*(h1 + h2 &
    + h3) + b21*b22*(2d0*h1*h2 + 6d0*h2**2d0 + h1*h3 &
    + 4d0*h2*h3 + 2d0*h3**2d0) + b22**2d0*(h1**2d0 &
    + 6d0*h2**2d0 + 3d0*h2*h3 + h3**2d0 + h1*(3d0*h2 + h3)))

  bs4t1(4,6)=(4d0/45d0)*((-b21**2d0)*(h1 - h2)*(h1 + h2 &
    + h3) - b22**2d0*(h1 - h3)*(h1 + h2 + h3) &
    + b21*b22*(4d0*h1**2d0 + 2d0*h2**2d0 + 3d0*h2*h3 &
    + 2d0*h3**2d0 + 2d0*h1*(h2 + h3)))

  bs4t1(5,1)=(1d0/45d0)*(b21 + b22)*(b22*(-3d0*h1**2d0 &
    + 3d0*h2**2d0 + 2d0*h2*h3 + h3**2d0 - h1*(2d0*h2 + h3)) &
    + b21*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 + 3d0*h3**2d0 &
    - h1*(h2 + 2d0*h3)))

  bs4t1(5,2)=(1d0/45d0)*b21*(b21*(-h1**2d0 + 3d0*h2**2d0 &
    + h1*(h2 - 2d0*h3) + 2d0*h2*h3 - 3d0*h3**2d0) &
    + b22*(h1**2d0 + 15d0*h2**2d0 + 6d0*h2*h3 + h3**2d0 &
    + h1*(6d0*h2 + h3)))

  bs4t1(5,3)=(-(1d0/45d0))*b22*(b22*(h1**2d0 + 2d0*h1*h2 &
    + 3d0*h2**2d0 - h1*h3 - 2d0*h2*h3 - 3d0*h3**2d0) &
    - b21*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs4t1(5,4)=(-(4d0/45d0))*((-b21**2d0)*(h1 - h2)*(h1 + h2 &
    + h3) + b21*b22*(2d0*h1*h2 + 6d0*h2**2d0 + h1*h3 &
    + 4d0*h2*h3 + 2d0*h3**2d0) + b22**2d0*(h1**2d0 &
    + 6d0*h2**2d0 + 3d0*h2*h3 + h3**2d0 + h1*(3d0*h2 + h3)))

  bs4t1(5,5)=(4d0/45d0)*(b21*b22*(h1**2d0 + 3d0*h2**2d0 &
    + 4d0*h2*h3 + 3d0*h3**2d0 + 2d0*h1*(h2 + h3)) &
    + b22**2d0*(h1**2d0 + 6d0*h2**2d0 + 3d0*h2*h3 + h3**2d0 &
    + h1*(3d0*h2 + h3)) + b21**2d0*(h1**2d0 + h2**2d0 &
    + 3d0*h2*h3 + 6d0*h3**2d0 + h1*(h2 + 3d0*h3)))

  bs4t1(5,6)=(-(4d0/45d0))*((-b22**2d0)*(h1 - h3)*(h1 + h2 &
    + h3) + b21*b22*(h1*h2 + 2d0*h2**2d0 + 2d0*h1*h3 &
    + 4d0*h2*h3 + 6d0*h3**2d0) + b21**2d0*(h1**2d0 + h2**2d0 &
    + 3d0*h2*h3 + 6d0*h3**2d0 + h1*(h2 + 3d0*h3)))

  bs4t1(6,1)=(-(1d0/45d0))*(b21 + b22)*(b22*(12d0*h1**2d0 &
    + 5d0*h1*h2 + 2d0*h2**2d0 + 4d0*h1*h3 + 3d0*h2*h3 &
    + 4d0*h3**2d0) + b21*(-3d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 &
    + 3d0*h3**2d0 - h1*(h2 + 2d0*h3)))

  bs4t1(6,2)=(-(1d0/45d0))*b21*(b22*(h1 - h3)*(2d0*h1 - h2 &
    + 2d0*h3) + b21*(-h1**2d0 + 3d0*h2**2d0 + h1*(h2 - 2d0*h3) &
    + 2d0*h2*h3 - 3d0*h3**2d0))

  bs4t1(6,3)=(-(1d0/45d0))*b22*(b22*(4d0*h1**2d0 + 3d0*h1*h2 &
    + 2d0*h2**2d0 + 4d0*h1*h3 + 5d0*h2*h3 + 12d0*h3**2d0) &
    + b21*(h1**2d0 + h2**2d0 + 6d0*h2*h3 + 15d0*h3**2d0 &
    + h1*(h2 + 6d0*h3)))

  bs4t1(6,4)=(4d0/45d0)*((-b21**2d0)*(h1 - h2)*(h1 + h2 &
    + h3) - b22**2d0*(h1 - h3)*(h1 + h2 + h3) &
    + b21*b22*(4d0*h1**2d0 + 2d0*h2**2d0 + 3d0*h2*h3 &
    + 2d0*h3**2d0 + 2d0*h1*(h2 + h3)))

  bs4t1(6,5)=(-(4d0/45d0))*((-b22**2d0)*(h1 - h3)*(h1 + h2 &
    + h3) + b21*b22*(h1*h2 + 2d0*h2**2d0 + 2d0*h1*h3 &
    + 4d0*h2*h3 + 6d0*h3**2d0) + b21**2d0*(h1**2d0 + h2**2d0 &
    + 3d0*h2*h3 + 6d0*h3**2d0 + h1*(h2 + 3d0*h3)))

  bs4t1(6,6)=(4d0/45d0)*(b21*b22*(-h1**2d0 + h2**2d0 &
    + 2d0*h1*h3 + 4d0*h2*h3 + 9d0*h3**2d0) &
    + b22**2d0*(4d0*h1**2d0 + h2**2d0 + 2d0*h2*h3 + 4d0*h3**2d0 &
    + 2d0*h1*(h2 + h3)) + b21**2d0*(h1**2d0 + h2**2d0 &
    + 3d0*h2*h3 + 6d0*h3**2d0 + h1*(h2 + 3d0*h3)))

end subroutine bsnq4Term1

subroutine bsnq4Term2(bs4t2,b11,b12,b21,b22,h1,h2,h3)
implicit none

  real(kind=8),intent(out)::bs4t2(6,6)
  real(kind=8),intent(in)::b11,b12,b21,b22,h1,h2,h3 

  bs4t2(1,1)=(1d0/120d0)*(b21 + b22)*(b21*(h1 - h2) &
    + b22*(h1 - h3))*(6d0*h1 + h2 + h3)

  bs4t2(1,2)=(1d0/360d0)*b21*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(6d0*h1 + 5d0*h2 + h3)

  bs4t2(1,3)=(1d0/360d0)*b22*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(6d0*h1 + h2 + 5d0*h3)

  bs4t2(1,4)=(-(1d0/90d0))*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b22*(2d0*h2 + h3) + b21*(6d0*h1 + 2d0*h2 + h3))

  bs4t2(1,5)=(1d0/90d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b22*(2d0*h2 + h3) + b21*(h2 + 2d0*h3))

  bs4t2(1,6)=(-(1d0/90d0))*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b21*(h2 + 2d0*h3) + b22*(6d0*h1 + h2 + 2d0*h3))

  bs4t2(2,1)=(1d0/360d0)*(b21 + b22)*(5d0*h1 + 6d0*h2 &
    + h3)*(b21*(-h1 + h2) + b22*(-h1 + h3))

  bs4t2(2,2)=(-(1d0/120d0))*b21*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(h1 + 6d0*h2 + h3)

  bs4t2(2,3)=(1d0/360d0)*b22*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(h1 + 6d0*h2 + 5d0*h3)

  bs4t2(2,4)=(1d0/90d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(6d0*b22*h2 + b21*(2d0*h1 + 6d0*h2 + h3))

  bs4t2(2,5)=(1d0/90d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(-6d0*b22*h2 + b21*(h1 + 2d0*h3))

  bs4t2(2,6)=(-(1d0/90d0))*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b22*(-h1 + h3) + b21*(h1 + 2d0*h3))

  bs4t2(3,1)=(1d0/360d0)*(b21 + b22)*(5d0*h1 + h2 &
    + 6d0*h3)*(b21*(-h1 + h2) + b22*(-h1 + h3))

  bs4t2(3,2)=(1d0/360d0)*b21*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(h1 + 5d0*h2 + 6d0*h3)

  bs4t2(3,3)=(1d0/120d0)*b22*(h1 + h2 + 6d0*h3)*(b21*(-h1 &
    + h2) + b22*(-h1 + h3))

  bs4t2(3,4)=(1d0/90d0)*(b21*(h1 - h2) - b22*(h1 &
    + 2d0*h2))*(b21*(h1 - h2) + b22*(h1 - h3))

  bs4t2(3,5)=(1d0/90d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b22*(h1 + 2d0*h2) - 6d0*b21*h3)

  bs4t2(3,6)=(1d0/90d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(6d0*b21*h3 + b22*(2d0*h1 + h2 + 6d0*h3))

  bs4t2(4,1)=(1d0/90d0)*(b21 + b22)*(b21*(h1 - h2) &
    + b22*(h1 - h3))*(6d0*h1 + 2d0*h2 + h3)

  bs4t2(4,2)=(-(1d0/90d0))*b21*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(2d0*h1 + 6d0*h2 + h3)

  bs4t2(4,3)=(1d0/90d0)*b22*(-2d0*h1 - 2d0*h2 + h3)*(b21*(-h1 &
    + h2) + b22*(-h1 + h3))

  bs4t2(4,4)=(2d0/45d0)*(b21*(-h1 + h2) + b22*(-h1 &
    + h3))*(b21*(h1 - h2) - b22*(2d0*h1 + 3d0*h2 + h3))

  bs4t2(4,5)=(2d0/45d0)*(b21*(-h1 + h2) + b22*(-h1 &
    + h3))*(b21*(h1 + h2 + h3) + b22*(2d0*h1 + 3d0*h2 + h3))

  bs4t2(4,6)=(2d0/45d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*((-b22)*(2d0*h1 + h2) + b21*(h1 + h2 + h3))

  bs4t2(5,1)=(1d0/90d0)*(b21 + b22)*(b21*(h1 - h2) &
    + b22*(h1 - h3))*(h1 - 2d0*(h2 + h3))

  bs4t2(5,2)=(-(1d0/90d0))*b21*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(h1 + 6d0*h2 + 2d0*h3)

  bs4t2(5,3)=(1d0/90d0)*b22*(h1 + 2d0*h2 + 6d0*h3)*(b21*(-h1 &
    + h2) + b22*(-h1 + h3))

  bs4t2(5,4)=(2d0/45d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b21*(2d0*h2 + h3) + b22*(h1 + 3d0*h2 + 2d0*h3))

  bs4t2(5,5)=(2d0/45d0)*(b21*(-h1 + h2) + b22*(-h1 &
    + h3))*(b22*(h1 + 3d0*h2 + 2d0*h3) + b21*(h1 &
    + 2d0*h2 + 3d0*h3))

  bs4t2(5,6)=(2d0/45d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b22*(h2 + 2d0*h3) + b21*(h1 + 2d0*h2 + 3d0*h3))

  bs4t2(6,1)=(1d0/90d0)*(b21 + b22)*(b21*(h1 - h2) &
    + b22*(h1 - h3))*(6d0*h1 + h2 + 2d0*h3)

  bs4t2(6,2)=(1d0/90d0)*b21*(-2d0*h1 + h2 - 2d0*h3)*(b21*(-h1 &
    + h2) + b22*(-h1 + h3))

  bs4t2(6,3)=(-(1d0/90d0))*b22*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(2d0*h1 + h2 + 6d0*h3)

  bs4t2(6,4)=(-(2d0/45d0))*(b21*(-h1 + h2) + b22*(-h1 &
    + h3))*((-b21)*(2d0*h1 + h3) + b22*(h1 + h2 + h3))

  bs4t2(6,5)=(-(2d0/45d0))*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b22*(h1 + h2 + h3) + b21*(2d0*h1 + h2 + 3d0*h3))

  bs4t2(6,6)=(2d0/45d0)*(b21*(h1 - h2) + b22*(h1 &
    - h3))*(b22*(-h1 + h3) + b21*(2d0*h1 + h2 + 3d0*h3))

end subroutine bsnq4Term2

subroutine cFluxMat(cxFlux,cyFlux,b11,b12,b21,b22)
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

end subroutine cFluxMat

subroutine remainMat(ax,ay,dMat,b11,b12,b21,b22,h1,h2,h3)
implicit none

  real(kind=8),intent(out)::ax(6,3),ay(6,3),dMat(3,3)
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


  !Taking aux var equation as
  !w = d/dy(h d(eta)/dx) + d/dx(h d(eta)/dy)
  ! dMat(1,1)=(1d0/3d0)*(b11 + b12)*(b21 + b22)*(h1 + h2 + h3)

  ! dMat(1,2)=(-(1d0/6d0))*(b12*b21 + b11*(2d0*b21 &
  !   + b22))*(h1 + h2 + h3)

  ! dMat(1,3)=(-(1d0/6d0))*(b11*b22 + b12*(b21 &
  !   + 2d0*b22))*(h1 + h2 + h3)

  ! dMat(2,1)=(-(1d0/6d0))*(b12*b21 + b11*(2d0*b21 &
  !   + b22))*(h1 + h2 + h3)

  ! dMat(2,2)=(1d0/3d0)*b11*b21*(h1 + h2 + h3)

  ! dMat(2,3)=(1d0/6d0)*(b12*b21 + b11*b22)*(h1 + h2 + h3)

  ! dMat(3,1)=(-(1d0/6d0))*(b11*b22 + b12*(b21 &
  !   + 2d0*b22))*(h1 + h2 + h3)

  ! dMat(3,2)=(1d0/6d0)*(b12*b21 + b11*b22)*(h1 + h2 + h3)

  ! dMat(3,3)=(1d0/3d0)*b12*b22*(h1 + h2 + h3)

  !Taking aux var equation as
  !w = d/dx(h d(eta)/dx) + d/dy(h d(eta)/dy)
  dMat(1,1)=(1d0/6d0)*((b11 + b12)**2d0 + (b21 &
    + b22)**2d0)*(h1 + h2 + h3)

  dMat(1,2)=(-(1d0/6d0))*(b11**2d0 + b11*b12 + b21*(b21 &
    + b22))*(h1 + h2 + h3)

  dMat(1,3)=(-(1d0/6d0))*(b11*b12 + b12**2d0 + b22*(b21 &
    + b22))*(h1 + h2 + h3)

  dMat(2,1)=(-(1d0/6d0))*(b11**2d0 + b11*b12 + b21*(b21 &
    + b22))*(h1 + h2 + h3)

  dMat(2,2)=(1d0/6d0)*(b11**2d0 + b21**2d0)*(h1 + h2 + h3)

  dMat(2,3)=(1d0/6d0)*(b11*b12 + b21*b22)*(h1 + h2 + h3)

  dMat(3,1)=(-(1d0/6d0))*(b11*b12 + b12**2d0 + b22*(b21 &
    + b22))*(h1 + h2 + h3)

  dMat(3,2)=(1d0/6d0)*(b11*b12 + b21*b22)*(h1 + h2 + h3)

  dMat(3,3)=(1d0/6d0)*(b12**2d0 + b22**2d0)*(h1 + h2 + h3)

end subroutine remainMat
