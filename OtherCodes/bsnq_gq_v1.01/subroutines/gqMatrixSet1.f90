subroutine gqMatrixSet1(npt,npl,nele,conn,ivl,ivq,linkl,&
  linkq,jacb,dep,por,stoneCnst,mass1,mass2,porLam,touMass,&
  gCxF,gCyF)
use shapeFnc
implicit none

  integer(kind=C_K1),intent(in)::npt,npl,nele,conn(nele,6)
  integer(kind=C_K1),intent(in)::ivl(0:npt),ivq(0:npt)
  integer(kind=C_K1),intent(in)::linkl(ivl(0)*npt)
  integer(kind=C_K1),intent(in)::linkq(ivq(0)*npt)
  integer(kind=C_K1)::gRow,gCol,lRow,lCol
  integer(kind=C_K1)::nlinkl(ivl(0)),nlinkq(ivq(0))

  type(jacbType)::jacb(nele)

  real(kind=C_K2),intent(in)::stoneCnst(10),por(npl)
  real(kind=C_K2),intent(in)::dep(npt)
  real(kind=C_K2),intent(out)::porLam(ivq(0)*npt)  
  real(kind=C_K2),intent(out)::touMass(ivl(0)*npt)
  real(kind=C_K2),intent(out)::mass1(ivq(0)*npt)
  real(kind=C_K2),intent(out)::mass2(ivl(0)*npl)  
  real(kind=C_K2),intent(out)::gCxF(ivq(0)*npl)
  real(kind=C_K2),intent(out)::gCyF(ivq(0)*npl)
  real(kind=C_K2)::lM1(6,6),lM2(3,3),lPLam(6,6)
  real(kind=C_K2)::lTouMass(6,3)
  real(kind=C_K2)::lCxF(3,6),lCyF(3,6)
  real(kind=C_K2)::jD(nGP),jD11(nGP),jD12(nGP)
  real(kind=C_K2)::jD21(nGP),jD22(nGP)
  real(kind=C_K2)::ePor(3),ed(3),eh(6)  
  real(kind=C_K2)::sPor,sD,sH
  real(kind=C_K2)::sc1

  ! lh = local depth still water
  ! ld = h+eta     

  mass1=0d0
  mass2=0d0
  porLam=0d0
  touMass=0d0 !ToBeDone
  gCxF=0d0  !Second term toBeDone
  gCyF=0d0  !Second term toBeDone
  
  do iel=1,nele
    nq=conn(iel,:)
    jD=jacb(iel)%D
    jD11=jacb(iel)%D11
    jD12=jacb(iel)%D12
    jD21=jacb(iel)%D21
    jD22=jacb(iel)%D22
    ePor=por(nq(1:3))
    eh=dep(nq)

    lM1=0d0    
    lM2=0d0
    lPLam=0d0
    lCxF=0d0
    lCyF=0d0

    do k=1,nGP
      !! Scalar      
      sPor=0d0      
      do i=1,3
        sPor=sPor+sh3F(i,k)*ePor(i)
      enddo
      sc1=(1d0-sPor)**3/sPor*stoneCnst(5)

      !! Element Matrices
      !6x6
      do i=1,6
        do j=1,6
          tmpr1=sh6F(i,k)*sh6F(j,k)*jD(k)
          tmpr2=tmpr1*sc1
          lM1(i,j)=lM1(i,j)+(gpW(k)*tmpr1)
          lPLam(i,j)=lPLam(i,j)+(gpW(k)*tmpr2)
        enddo
      enddo

      !3x3
      do i=1,3
        do j=1,3
          tmpr1=sh3F(i,k)*sh3F(j,k)*jD(k)
          lM2(i,j)=lM2(i,j)+(gpW(k)*tmpr1)
        enddo
      enddo

      !3x6
      do i=1,3
        do j=1,6          
          tmpr1=-sh3F(i,k)*(sh6FE(j,k)*jD11(k) &
            + sh6FN(j,k)*jD12(k))*jD(k)/sPor
          lCxF(i,j)=lCxF(i,j)+(gpW(k)*tmpr1)

          tmpr1=-sh3F(i,k)*(sh6FE(j,k)*jD21(k) &
            + sh6FN(j,k)*jD22(k))*jD(k)/sPor
          lCyF(i,j)=lCyF(i,j)+(gpW(k)*tmpr1)
        enddo
      enddo

      !! Storages
      !6x6
      do lRow=1,6
        gRow=nq(lRow)
        i=(gRow-1)*ivq(0)
        nlinkq=linkq(i+1:i+ivq(0))
        do lCol=1,6
          gCol=nq(lCol)
          do j=1,ivq(gRow)
            if(nlinkq(j).eq.gCol) exit
          enddo
          j=i+j
          mass1(j)=mass1(j)+lM1(lRow,lCol)
          porLam(j)=porLam(j)+lPLam(lRow,lCol)
        enddo
      enddo

      !3x3
      do lRow=1,3
        gRow=nq(lRow)
        i=(gRow-1)*ivl(0)
        nlinkl=linkl(i+1:i+ivl(0))
        do lCol=1,3
          gCol=nq(lCol)
          do j=1,ivl(gRow)
            if(nlinkl(j).eq.gCol) exit
          enddo
          j=i+j
          mass2(j)=mass2(j)+lM2(lRow,lCol)          
        enddo
      enddo

      !3x6
      do lRow=1,3
        gRow=nq(lRow)
        i=(gRow-1)*ivq(0)
        nlinkq=linkq(i+1:i+ivq(0))
        do lCol=1,6
          gCol=nq(lCol)
          do j=1,ivq(gRow)
            if(nlinkq(j).eq.gCol) exit
          enddo
          j=i+j
          gCxF(j)=gCxF(j)+lCxF(lRow,lCol)
          gCyF(j)=gCyF(j)+lCyF(lRow,lCol)
        enddo
      enddo

    enddo
  enddo

end subroutine gqMatrixSet1