!!-----------------------------setMFree----------------------------!!
  subroutine setMFree(b)
  implicit none

    class(bsnqCase),intent(inout)::b
    integer(kind=C_K1)::nn,err
    integer(kind=C_K1),allocatable::neid(:),newrk(:)
    real(kind=C_K2)::cx,cy,rad
    real(kind=C_K2),allocatable::nedr(:),phi(:),phiDx(:),phiDy(:)

    call system_clock(b%sysC(5))
    allocate(b%mfl(b%npl))
    allocate(neid(b%npt), newrk(b%npt), nedr(b%npt))
    allocate(phi(b%npt),phiDx(b%npt),phiDy(b%npt))

    tmpr5=1.2d0 !Coeff of multipliciation to max rad in linkList
    do i=1,b%npl

      call findRadLinkList(i, b%npt, b%Sz(4), b%ivq, b%linkq, b%cor, &
        tmpr5, rad, i2, j, k1, tmpr3, cx, cy)

      call findNeiLinkList(i, rad, b%npt, b%Sz(4), b%ivq, &
        b%linkq, b%cor, b%npt, nn, neid, newrk, nedr)

      call mls2DDx(cx, cy, nn, rad, b%cor(neid(1:nn),1), &
        b%cor(neid(1:nn),2), phi(1:nn), phiDx(1:nn), phiDy(1:nn), err)

      if(err.ne.0)then
        write(mf,'(" [ERR] No MFree at node ", I10)')i
        write(mf,'(" [---] Cx, Cy ",2F15.6)')cx,cy
      endif

      call b%mfl(i)%setPoi(nn, nn, i, cx, cy, rad, neid(1:nn), &
        phi(1:nn), phiDx(1:nn), phiDy(1:nn))

      ! if(i.ne.1040) cycle
      ! write(*,'(I15,2F15.6)')i,b%mfl(i)%cx,b%mfl(i)%cy
      ! write(*,'(2I15,F15.6)')b%mfl(i)%nn,b%mfl(i)%nnMax,b%mfl(i)%rad
      ! do j=1,nn
      !   write(*,'(I15,3F15.6)')b%mfl(i)%neid(j), b%mfl(i)%phi(j), &
      !     b%mfl(i)%phiDx(j), b%mfl(i)%phiDy(j)
      ! enddo
      ! write(*,'(F15.6)')sum(b%mfl(i)%phi(1:nn))

    enddo  

    deallocate(neid,newrk,nedr,phi,phiDx,phiDy)

    call system_clock(b%sysC(6))
    write(9,*)"[MSG] Done setMFree"
    write(9,'(" [TIM] ",F15.4)')1d0*(b%sysC(6)-b%sysC(5))/b%sysRate
    write(9,*)    

  end subroutine setMFree
!!---------------------------End setMFree--------------------------!!
