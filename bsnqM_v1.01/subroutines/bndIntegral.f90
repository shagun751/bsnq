subroutine bndIntegral1(npl,npt,nele,nbnd,conn,mabnd,ivl,ivq,&
  linkl,linkq,invJ,bndS,dep,gFp,gFq)
use bsnqGlobVars
implicit none

  integer(kind=C_K1),intent(in)::npl,npt,nele,nbnd,conn(nele,6)
  integer(kind=C_K1),intent(in)::ivl(0:npt),ivq(0:npt)
  integer(kind=C_K1),intent(in)::linkl(ivl(0)*npt)
  integer(kind=C_K1),intent(in)::linkq(ivq(0)*npt)
  integer(kind=C_K1),intent(in)::mabnd(nbnd,6)
  integer(kind=C_K1)::i,j,k,ele,en(6),gRow,gCol,lRow,lCol
  integer(kind=C_K1)::nlinkl(ivl(0)),nlinkq(ivq(0))
  
  real(kind=C_K2),intent(in)::invJ(nele,5),bndS(nbnd,3),dep(npt)
  real(kind=C_K2),intent(out)::gFp(ivq(0)*npt),gFq(ivq(0)*npt)
  real(kind=C_K2)::cnst(10)
  real(kind=C_K2)::nx,ny,sL
  real(kind=C_K2)::lFp(6,6),lFq(6,6)
  real(kind=C_K2)::b11,b12,b21,b22


  cnst(1)=grav                    !gravity
  cnst(2)=BsqC                    !Bsnq constant
  cnst(3)=BsqC+(1d0/3d0)          !dervied 1
  cnst(4)=2d0*BsqC+(1d0/3d0)      !derived 2
  cnst(5)=2d0*BsqC+(1d0/2d0)      !derived 3
  cnst(6)=grav*BsqC               !derived 4

  gFp=0d0
  gFq=0d0

  do i=1,nbnd
    ele=mabnd(i,3)
    en=conn(ele,:)

    b11=invJ(ele,1)
    b12=invJ(ele,2)
    b21=invJ(ele,3)
    b22=invJ(ele,4)
    nx=bndS(i,1)
    ny=bndS(i,2)
    sL=bndS(i,3)

    if(mabnd(i,6).eq.1) then
      call bsnqBnd12(lFp,lFq,cnst,b11,b12,b21,b22,nx,ny,&
        dep(en(1)),dep(en(2)),dep(en(3)))
    else if(mabnd(i,6).eq.2) then
      call boundaryCase23(lFp,lFq,lFn,lFw,nx,ny,lp,lq,lw,&
        leta,lh,cnst,b11,b12,b21,b22)
    else if(mabnd(i,6).eq.3) then
      call boundaryCase31(lFp,lFq,lFn,lFw,nx,ny,lp,lq,lw,&
        leta,lh,cnst,b11,b12,b21,b22)

  enddo

end subroutine bndIntegral1



subroutine bsnqBnd12(lFp,lFq,cnst,b11,b12,b21,b22,nx,ny,&
  h1,h2,h3)
use bsnqGlobVars
implicit none
  
  real(kind=C_K2),intent(out)::lFp(6,6),lFq(6,6)
  real(kind=C_K2),intent(in)::cnst(10),b11,b12,b21,b22
  real(kind=C_K2),intent(in)::nx,ny,h1,h2,h3
  real(kind=C_K2)::tmp(6,6)

  call bsnqTermBndInt(tmp,b11,b12,b21,b22)

end subroutine