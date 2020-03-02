!!-----------------------------setMFree----------------------------!!
  subroutine setMFree(b)
  implicit none

    class(bsnqCase),intent(inout)::b
    integer(kind=C_K1)::nn,err
    integer(kind=C_K1),allocatable::neid(:),newrk(:)
    real(kind=C_K2)::cx,cy,rad
    real(kind=C_K2),allocatable::nedr(:),phi(:),phiDx(:),phiDy(:)

    call system_clock(b%sysC(5))
    allocate(b%pObf(b%npt))
    allocate(neid(b%npt), newrk(b%npt), nedr(b%npt))
    allocate(phi(b%npt),phiDx(b%npt),phiDy(b%npt))

    tmpr5=1.2d0 !Coeff of multipliciation to max rad in linkList
    do i=1,b%npt

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

      call b%pObf(i)%setPoi(nn, nn, i, cx, cy, rad, neid(1:nn), &
        phi(1:nn), phiDx(1:nn), phiDy(1:nn))

      ! if(i.ne.1040) cycle
      ! write(*,'(I15,2F15.6)')i,b%pObf(i)%cx,b%pObf(i)%cy
      ! write(*,'(2I15,F15.6)')b%pObf(i)%nn,b%pObf(i)%nnMax,b%pObf(i)%rad
      ! do j=1,nn
      !   write(*,'(I15,3F15.6)')b%pObf(i)%neid(j), b%pObf(i)%phi(j), &
      !     b%pObf(i)%phiDx(j), b%pObf(i)%phiDy(j)
      ! enddo
      ! write(*,'(F15.6)')sum(b%pObf(i)%phi(1:nn))

    enddo  

    deallocate(neid,newrk,nedr,phi,phiDx,phiDy)

    call b%bDf%init( b%npt )

    call system_clock(b%sysC(6))
    write(9,*)"[MSG] Done setMFree"
    write(9,'(" [TIM] ",F15.4)')1d0*(b%sysC(6)-b%sysC(5))/b%sysRate
    write(9,*)    

  end subroutine setMFree
!!---------------------------End setMFree--------------------------!!



!!-----------------------------calcDerv----------------------------!!
  subroutine calcDerv(b)
  implicit none

    class(bsnqCase),intent(inout)::b

    b%ur = b%tOb(0)%p / b%tOb(0)%tD
    b%vr = b%tOb(0)%q / b%tOb(0)%tD

    !! First derivative
    !$OMP PARALLEL DEFAULT(shared) &
    !$OMP   PRIVATE(i,i2,j,k2)
    !$OMP DO SCHEDULE(dynamic,100)    
    do i=1,b%npt
      i2 = b%pObf(i)%bsnqId
      k2 = b%pObf(i)%nn

      call calcGrad( k2, b%pObf(i)%phiDx, &
        b%ur( b%pObf(i)%neid ), b%bDf%ux( i2 ), j )

      call calcGrad( k2, b%pObf(i)%phiDx, &
        b%tOb(0)%p( b%pObf(i)%neid ), b%bDf%px( i2 ), j )

    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    !! Second derivative
    !$OMP PARALLEL DEFAULT(shared) &
    !$OMP   PRIVATE(i,i2,j,k2)
    !$OMP DO SCHEDULE(dynamic,100)    
    do i=1,b%npt
      i2 = b%pObf(i)%bsnqId
      k2 = b%pObf(i)%nn

      call calcGrad( k2, b%pObf(i)%phiDx, &
        b%bDf%ux( b%pObf(i)%neid ), b%bDf%uxx( i2 ), j )

      call calcGrad( k2, b%pObf(i)%phiDx, &
        b%bDf%px( b%pObf(i)%neid ), b%bDf%pxx( i2 ), j )

    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    !! Third derivative
    !$OMP PARALLEL DEFAULT(shared) &
    !$OMP   PRIVATE(i,i2,j,k2)
    !$OMP DO SCHEDULE(dynamic,100)    
    do i=1,b%npt
      i2 = b%pObf(i)%bsnqId
      k2 = b%pObf(i)%nn

      call calcGrad( k2, b%pObf(i)%phiDx, &
        b%bDf%uxx( b%pObf(i)%neid ), b%bDf%uxxx( i2 ), j )

      call calcGrad( k2, b%pObf(i)%phiDx, &
        b%bDf%pxx( b%pObf(i)%neid ), b%bDf%pxxx( i2 ), j )

    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

  end subroutine calcDerv
!!---------------------------End calcDerv--------------------------!!
