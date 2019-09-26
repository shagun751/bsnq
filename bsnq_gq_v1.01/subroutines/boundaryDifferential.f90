subroutine etaBC(npoinl,npoint,nelem,nbndpoi,maxNeEle,conn,&
	poi2ele,bndNode,ivl,linkl,invJ,bndNormal,gdNdn)
implicit none

  integer(kind=4),intent(in)::npoinl,npoint,nelem,nbndpoi
  integer(kind=4),intent(in)::maxNeEle,conn(nelem,6)
  integer(kind=4),intent(in)::poi2ele(npoint,maxNeEle)
  integer(kind=4),intent(in)::bndNode(nbndpoi),ivl(0:npoint)
  integer(kind=4),intent(in)::linkl(ivl(0)*npoint)
  integer(kind=4)::i,j,k,i2,j2,ele,n(3),poi,gCol,nlinkl(ivl(0))

  real(kind=8),intent(in)::invJ(nelem,5),bndNormal(npoint,2)
  real(kind=8),intent(out)::gdNdn(ivl(0)*npoinl)
  real(kind=8)::ldNdx(3),ldNdy(3),ldNdn(3),nx,ny

  gdNdn=0

  do j2=1,nbndpoi
    poi=bndNode(j2)
    if(poi.gt.npoinl) goto 13
    do i=1,maxNeEle
      ele=poi2ele(poi,i)
      if(ele.eq.0) goto 11
      n=conn(ele,1:3)

      ldNdx(1)=-(invJ(ele,1)+invJ(ele,2))
      ldNdx(2)=invJ(ele,1)
      ldNdx(3)=invJ(ele,2)

      ldNdy(1)=-(invJ(ele,3)+invJ(ele,4))
      ldNdy(2)=invJ(ele,3)
      ldNdy(3)=invJ(ele,4)

      nx=bndNormal(poi,1)
      ny=bndNormal(poi,2)

      ldNdn=((ldNdx*nx)+(ldNdy*ny))*invJ(ele,5)/6.0d0

      k=(poi-1)*ivl(0)
      nlinkl=linkl(k+1:k+ivl(0))
      do j=1,3
        gCol=n(j)
        do i2=1,ivl(poi)
          if(nlinkl(i2).eq.gCol) goto 12
        enddo
        write(9,*)"[Err] node conn missing in dNdn at",poi
        stop
        12 gdNdn(k+i2)=gdNdn(k+i2)+ldNdn(j)
      enddo
    enddo
    11 continue
  enddo
  13 continue

  ! do i=1,npoinl
  !   k=(i-1)*ivl(0)
  !   write(9,*)i,":",gdNdn(k+1:k+ivl(i))    
  ! enddo

end subroutine

! subroutine pqBC(npoinl,npoint,nelem,nbndpoi,maxNeEle,conn,&
!   poi2ele,bndNode,ivq,linkq,invJ,bndNormal)