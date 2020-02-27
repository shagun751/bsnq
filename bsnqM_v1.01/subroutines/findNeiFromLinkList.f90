!!-------------------------findRadLinkList-------------------------!!
  subroutine findRadLinkList(ip,npt,szJv,ivq,jvq,cor,coef,&
    rad,i2,j,k1,dr,cx,cy)
  use bsnqGlobVars
  implicit none

    integer(kind=C_K1),intent(in)::ip,npt,szJv
    integer(kind=C_K1),intent(in)::ivq(0:npt),jvq(szJv)
    real(kind=C_K2),intent(in)::cor(npt,2),coef

    integer(kind=C_K1),intent(out)::j,k1,i2
    real(kind=C_K2),intent(out)::rad,dr,cx,cy

    ! Radius
    rad=0d0
    cx=cor(ip,1)
    cy=cor(ip,2)
    k1=(ip-1)*ivq(0)
    do j=k1+1,k1+ivq(ip)
      i2=jvq(j)      
      dr=(cor(i2,1)-cx)**2 + (cor(i2,2)-cy)**2
      if(rad.lt.dr) rad=dr      
    enddo
    if(rad.lt.1e-10)then
      write(9,'(" [ERR] Check radius calculation at node",I10)')ip
      stop
    endif    
    rad=sqrt(rad)*coef

  end subroutine findRadLinkList
!!-----------------------End findRadLinkList-----------------------!!




!!-------------------------findNeiLinkList-------------------------!!
  subroutine findNeiLinkList(ip,rad,npt,szJv,ivq,jvq,cor,&
    nnMax,nn,neid,newrk,nedr)
  use bsnqGlobVars
  implicit none

    !! find nearest neighbours in radius 'rad' using FEM link list

    integer(kind=C_K1),intent(in)::ip,npt,szJv,nnMax
    integer(kind=C_K1),intent(in)::ivq(0:npt),jvq(szJv)
    real(kind=C_K2),intent(in)::rad,cor(npt,2)

    integer(kind=C_K1),intent(out)::nn,neid(nnMax),newrk(nnMax)
    real(kind=C_K2),intent(out)::nedr(nnMax)

    integer(kind=C_K1)::i2,i3,j,j2,j3,k1,k2,k3
    real(kind=C_K2)::px,py,rad2,rat2

    nn=0

    px=cor(ip,1)
    py=cor(ip,2)
    rad2=rad**2

    k2=0
    k3=1
    newrk(1)=ip
    nedr(1)=0d0

    301 continue
    k1=k2+1
    k2=k3
    do j=k1,k2
      i2=newrk(j)
      i3=(i2-1)*ivq(0)
      do j2=1,ivq(i2)
        if(count(newrk(1:k3).eq.jvq(i3+j2)).eq.0)then
          j3=jvq(i3+j2)
          rat2 = ((cor(j3,1)-px)**2 + (cor(j3,2)-py)**2)/rad2
          k3=k3+1
          newrk(k3)=j3
          nedr(k3)=rat2
        endif
      enddo
    enddo

    if(maxval(nedr(k2+1:k3)).lt.1d0) goto 301 !! looping till 2x radius

    do i2=1,k3
      if(nedr(i2).gt.1d0)cycle
      nn=nn+1
      neid(nn)=newrk(i2)
    enddo

    ! write(*,'(I10)')nn
    ! write(*,'(2F10.4)')cor(ip,:)

    ! do i2=1,nn
    !   j3=neid(i2)
    !   rat2 = sqrt(((cor(j3,1)-px)**2 + (cor(j3,2)-py)**2)/rad2)
    !   write(*,'(I10,3F10.4)')j3,rat2,cor(j3,:)
    ! enddo


  end subroutine findNeiLinkList
!!-----------------------End findNeiLinkList-----------------------!!



