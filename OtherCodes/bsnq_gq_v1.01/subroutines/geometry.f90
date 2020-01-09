subroutine jacbInvLin(npoint,nelem,conn,coord,invJ)
implicit none
  
  ![Note]:
  !invJ(i,5)=invJ11,invJ12,invJ21,invJ22,|J|

  integer(kind=4),intent(in)::nelem,npoint
  integer(kind=4),intent(in)::conn(nelem,6)
  integer(kind=4)::i,j,k,n(3)

  real(kind=8),intent(in)::coord(npoint,2)
  real(kind=8),intent(out)::invJ(nelem,5)
  real(kind=8)::temp(5)

  invJ=0

  do i=1,nelem
    n=(/ (conn(i,j),j=1,3) /)

    temp(1)=coord(n(3),2)-coord(n(1),2)
    temp(2)=coord(n(1),2)-coord(n(2),2)
    temp(3)=coord(n(1),1)-coord(n(3),1)
    temp(4)=coord(n(2),1)-coord(n(1),1)
    temp(5)=(temp(1)*temp(4)-temp(2)*temp(3))

    if(temp(5).le.0) then
      write(9,*)"[Err] Invld Ele!! Clk nodes or 0 area at ele",i
      stop
    endif
    
    temp(1:4)=temp(1:4)/temp(5)    

    invJ(i,:)=temp
  enddo

end subroutine jacbInvLin


subroutine jacbInvQuad(npt,nele,conn,coord,jacb)
use shapeFnc
implicit none

  integer(kind=C_K1),intent(in)::npt,nele,conn(nele,6)
  real(kind=C_K2),intent(in)::coord(npt,2)  
  real(kind=C_K2)::ecx(6),ecy(6)
  type(jacbType),intent(out)::jacb(nele)

  !! Note:
  ! jacb(i)%D(j) = value of D inside element i at integration point j

  do iel=1,nele
    nq=conn(iel,:)    
    ecx=coord(nq,1)
    ecy=coord(nq,2)

    do i=1,nGP !Looping over integration points
      tmpr1=0d0
      tmpr2=0d0
      tmpr3=0d0
      tmpr4=0d0
      do j=1,6 !Looping over element nodes
        tmpr1=tmpr1+sh6FE(j,i)*ecx(j)
        tmpr2=tmpr2+sh6FE(j,i)*ecy(j)
        tmpr3=tmpr3+sh6FN(j,i)*ecx(j)
        tmpr4=tmpr4+sh6FN(j,i)*ecy(j)
      enddo
      tmpr5=dabs(tmpr1*tmpr4 - tmpr2*tmpr3)
      if(tmpr5.lt.1d-10)then
        write(9,*)"[ERR] Jacobian determinant < 1e-10, ele",iel
        write(9,*)"[---] Possible Clk nodes"
        stop
      endif
      jacb(iel)%D11(i)=tmpr4/tmpr5
      jacb(iel)%D12(i)=-tmpr2/tmpr5
      jacb(iel)%D21(i)=-tmpr3/tmpr5
      jacb(iel)%D22(i)=tmpr1/tmpr5            
      jacb(iel)%D(i)=tmpr5
    enddo
  enddo

  ! !! Debug comments
  ! write(tf,*)"[DBG] Jacobian determinant"
  ! do iel=1,nele
  !   write(tf,'(I6,9F15.6)')iel,jacb(iel)%D(:)
  !   tmpr1=0d0
  !   do i=1,9
  !     tmpr1=tmpr1+jacb(iel)%D(i)*shW(i)      
  !   enddo
  !   write(tf,*)iel,tmpr1
  ! enddo  
end subroutine jacbInvQuad


subroutine bndSideInfo(npoinl,npoint,nelem,nbnd,coord,mabnd,&
  bndSide)
implicit none
  
  integer(kind=4),intent(in)::npoinl,npoint,nelem,nbnd
  integer(kind=4),intent(in)::mabnd(nbnd,6)
  integer(kind=4)::i,j,k,n1,n2

  real(kind=8),intent(in)::coord(npoint,2)
  real(kind=8),intent(out)::bndSide(nbnd,3)
  real(kind=8)::x1,x2,y1,y2,s,nx,ny
  
  do i=1,nbnd
    n1=mabnd(i,1)
    n2=mabnd(i,2)
    x1=coord(n1,1)
    x2=coord(n2,1)
    y1=coord(n1,2)
    y2=coord(n2,2)

    nx=y2-y1
    ny=x1-x2
    s=dsqrt(nx**2 + ny**2)
    nx=nx/s
    ny=ny/s

    bndSide(i,1)=nx
    bndSide(i,2)=ny
    bndSide(i,3)=s
  enddo
end subroutine

subroutine bndNodeNormal(npoinl,npoint,nelem,nbnd,coord,mabnd,&
  bndSide,bndNormal)
implicit none
  
  integer(kind=4),intent(in)::npoint,npoinl,nelem,nbnd
  integer(kind=4),intent(in)::mabnd(nbnd,6)
  integer(kind=4)::i,j,k

  real(kind=8),intent(in)::coord(npoint,2),bndSide(nbnd,3)
  real(kind=8),intent(out)::bndNormal(npoint,2)
  real(kind=8)::temp

  bndNormal=0

  do i=1,nbnd
    bndNormal(mabnd(i,1),1)=bndNormal(mabnd(i,1),1)+&
      (bndSide(i,1)*bndSide(i,3))
    bndNormal(mabnd(i,2),1)=bndNormal(mabnd(i,2),1)+&
      (bndSide(i,1)*bndSide(i,3))
    bndNormal(mabnd(i,5),1)=bndNormal(mabnd(i,5),1)+&
      (bndSide(i,1)*bndSide(i,3))

    bndNormal(mabnd(i,1),2)=bndNormal(mabnd(i,1),2)+&
      (bndSide(i,2)*bndSide(i,3))
    bndNormal(mabnd(i,2),2)=bndNormal(mabnd(i,2),2)+&
      (bndSide(i,2)*bndSide(i,3))
    bndNormal(mabnd(i,5),2)=bndNormal(mabnd(i,5),2)+&
      (bndSide(i,2)*bndSide(i,3))
  enddo

  do i=1,npoint
    temp=dsqrt(bndNormal(i,1)**2 + bndNormal(i,2)**2)
    if(temp.ne.0) bndNormal(i,:)=bndNormal(i,:)/temp
    !write(9,*)i,":",bndNormal(i,1),bndNormal(i,2)
  enddo

end subroutine

subroutine middleNode(npoinl,npoinq,npoint,nelem,nedge,maxNePoi,&
    coord,depth,conn,poi2poi,npoisur)
implicit none

  integer(kind=4),intent(inout)::npoinl,npoinq,npoint,nelem
  integer(kind=4),intent(in)::maxNePoi,nedge
  integer(kind=4),intent(inout)::conn(nelem,6),npoisur(npoint,3)
  integer(kind=4),intent(inout)::poi2poi(npoint,maxNePoi)
  integer(kind=4)::i,j,k,l,ln1,ln2,i2,j2,k2,qp

  real(kind=8),intent(inout)::coord(npoint,2),depth(npoint)
  real(kind=8)::qx,qy,qd

  npoinq=0
  npoint=npoinl

  do i=1,nelem
    do j=1,3
      ln1=conn(i,j)
      ln2=conn(i,mod(j,3)+1)

      qp=0
      !Search for common middle node
      do i2=1,npoisur(ln1,1)
        do j2=1,npoisur(ln2,1)
          if(poi2poi(ln1,i2).eq.poi2poi(ln2,j2)) then
            goto 11
          endif
        enddo
      enddo

      !if common middle node not found
      npoint=npoint+1
      if(npoint.gt.(nedge+npoinl)) then
        write(9,*)"[Err] npoint .NE. npoinl + nedge"
        stop
      endif
      qp=npoint
      coord(qp,1)=0.5d0*(coord(ln1,1)+coord(ln2,1))
      coord(qp,2)=0.5d0*(coord(ln1,2)+coord(ln2,2))
      depth(qp)=0.5d0*(depth(ln1)+depth(ln2))
      npoisur(ln1,1)=npoisur(ln1,1)+1
      poi2poi(ln1,npoisur(ln1,1))=qp
      npoisur(ln2,1)=npoisur(ln2,1)+1
      poi2poi(ln2,npoisur(ln2,1))=qp
      goto 12

      !if common middle node is found
      11 qp=poi2poi(ln1,i2)

      12 conn(i,j+3)=qp

    enddo
  enddo

  npoinq=npoint-npoinl

end subroutine middleNode

subroutine bndMiddleNodes(npoinl,npoint,nelem,nbnd,conn,mabnd)
implicit none
  
  integer(kind=4),intent(in)::npoinl,npoint,nelem,nbnd
  integer(kind=4),intent(in)::conn(nelem,6)
  integer(kind=4),intent(inout)::mabnd(nbnd,6)
  integer(kind=4)::i,j,k,ele,p1,p2,p3,n(3)
  integer(kind=4)::temp(3),temp2(3)

  temp=(/ 4,6,5 /)
  temp2=(/ 2,3,1 /)

  do i=1,nbnd
    n(1)=mabnd(i,1)
    n(2)=mabnd(i,2)
    ele=mabnd(i,3)    

    do j=1,3
      if(conn(ele,j).eq.n(1)) goto 11
    enddo
    write(9,*)"[Err] Node not found at boundary for ele", ele
    stop
    11 p1=j

    do j=1,3
      if(conn(ele,j).eq.n(2)) goto 12
    enddo
    write(9,*)"[Err] Node not found at boundary for ele", ele
    stop
    12 p2=j

    if(p1.eq.p2) then
      write(9,*)"[Err] Both boundary nodes are same for ele", ele
      stop
    endif

    if(p2.ne.temp2(p1)) then
      write(9,*)"[Err] Side node clk for ele", ele
      stop
    endif

    p3=temp(p1+p2-2)
    n(3)=conn(ele,p3)
    mabnd(i,5)=n(3)
    mabnd(i,6)=p1    

  enddo
end subroutine bndMiddleNodes
