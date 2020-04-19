!!--------------------------findEleForLocXY------------------------!!
  subroutine findEleForLocXY1(b,xin,yin,eleOut,ep,et)
  implicit none

    class(bsnqCase),intent(in)::b
    real(kind=C_K2),intent(in)::xin,yin
    integer(kind=C_K1),intent(out)::eleOut
    real(kind=C_K2),intent(out)::ep,et
    integer(kind=C_K1)::ind(3,2),p1,p2,dir
    real(kind=C_K2)::xy(3,2),vec1(2),vec2(2),res

    ind(1,:)=(/ 1,2 /)
    ind(2,:)=(/ 2,3 /)
    ind(3,:)=(/ 3,1 /)

    do iel = 1, b%nele
      nl=b%conn(iel,1:3)
      do i=1,3
        xy(i,:)=b%cor(nl(i),:)
      enddo

      dir=1
      do i=1,3
        p1=ind(i,1)
        p2=ind(i,2)
        vec1(1)=xy(p2,1)-xy(p1,1)
        vec1(2)=xy(p2,2)-xy(p1,2)
        vec2(1)=xin-xy(p1,1)
        vec2(2)=yin-xy(p1,2)
        call vecCross2D(vec1(1),vec1(2),vec2(1),vec2(2),res)
        if(res.lt.0d0) dir=-1        
      enddo

      if(dir.eq.1)then
        eleOut=iel
        ep = ( b%invJ(iel,1)*(xin-xy(1,1)) ) &
          + ( b%invJ(iel,3)*(yin-xy(1,2)) )
        et = ( b%invJ(iel,2)*(xin-xy(1,1)) ) &
          + ( b%invJ(iel,4)*(yin-xy(1,2)) )
        return
      endif

    enddo

    eleOut=-1
    ep=-1
    et=-1

  end subroutine findEleForLocXY1

  subroutine findEleForLocXY2(b,np,xin,yin,eleOut,ep,et)
  implicit none
    
    class(bsnqCase),intent(in)::b
    integer(kind=C_K1),intent(in)::np
    real(kind=C_K2),intent(in)::xin(np),yin(np)
    integer(kind=C_K1),intent(out)::eleOut(np)
    real(kind=C_K2),intent(out)::ep(np),et(np)
    integer(kind=C_K1)::ind(3,2),p1,p2,dir,pp
    real(kind=C_K2)::xy(3,2),vec1(2),vec2(2),res,xp,yp

    ind(1,:)=(/ 1,2 /)
    ind(2,:)=(/ 2,3 /)
    ind(3,:)=(/ 3,1 /)


    !$OMP PARALLEL DEFAULT(shared) &
    !$OMP   PRIVATE(pp,xp,yp,iel,nl,xy,dir,i,p1,p2,vec1,vec2,res)
    !$OMP DO SCHEDULE(dynamic,10)
    do pp=1,np
      xp=xin(pp)
      yp=yin(pp)

      do iel = 1, b%nele
        nl=b%conn(iel,1:3)
        do i=1,3
          xy(i,:)=b%cor(nl(i),:)
        enddo

        dir=1
        do i=1,3
          p1=ind(i,1)
          p2=ind(i,2)
          vec1(1)=xy(p2,1)-xy(p1,1)
          vec1(2)=xy(p2,2)-xy(p1,2)
          vec2(1)=xp-xy(p1,1)
          vec2(2)=yp-xy(p1,2)
          call vecCross2D(vec1(1),vec1(2),vec2(1),vec2(2),res)
          if(res.lt.0d0) dir=-1        
        enddo

        if(dir.eq.1)then
          eleOut(pp)=iel
          ep(pp) = ( b%invJ(iel,1)*(xp-xy(1,1)) ) &
            + ( b%invJ(iel,3)*(yp-xy(1,2)) )
          et(pp) = ( b%invJ(iel,2)*(xp-xy(1,1)) ) &
            + ( b%invJ(iel,4)*(yp-xy(1,2)) )
          exit
        endif
      enddo
      
      if(dir.eq.1) cycle

      eleOut(pp)=-1
      ep(pp)=-1
      et(pp)=-1
    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL


  end subroutine findEleForLocXY2
!!------------------------End findEleForLocXY----------------------!!



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
    b%uhr = b%ur * b%dep
    b%vhr = b%vr * b%dep
    b%bDf%u = b%ur

    ! Storing ! d(U)/dx and d(Uh)/dx at t(n-1), t(n-2)    
    b%bDf%uxtn(:,2) = b%bDf%uxtn(:,1)  
    b%bDf%uhxtn(:,2) = b%bDf%uhxtn(:,1)    
    b%bDf%uxtn(:,1) = b%bDf%ux
    b%bDf%uhxtn(:,1) = b%bDf%uhx

    !b%ur = dsin(b%cor(:,1))

    !! First derivative
    !$OMP PARALLEL DEFAULT(shared) &
    !$OMP   PRIVATE(i,i2,j,k,k2)
    !$OMP DO SCHEDULE(dynamic,100)    
    do i=1,b%npt
      i2 = b%pObf(i)%bsnqId
      k2 = b%pObf(i)%nn

      call calcGrad( k2, b%pObf(i)%phiDx, &
        b%ur( b%pObf(i)%neid ), b%bDf%ux( i2 ), j )

      call calcGrad( k2, b%pObf(i)%phiDx, &
        b%uhr( b%pObf(i)%neid ), b%bDf%uhx( i2 ), j )

      call calcGrad( k2, b%pObf(i)%phiDx, &
        b%dep( b%pObf(i)%neid ), b%bDf%hx( i2 ), j )

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
        b%bDf%uhx( b%pObf(i)%neid ), b%bDf%uhxx( i2 ), j )

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
        b%bDf%uhxx( b%pObf(i)%neid ), b%bDf%uhxxx( i2 ), j )

    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL    

    b%bDf%uxt = ( 3d0*b%bDf%ux - 4d0*b%bDf%uxtn(:,1) &
      + b%bDf%uxtn(:,2) )/2d0/b%dt
    b%bDf%uhxt = ( 3d0*b%bDf%uhx - 4d0*b%bDf%uhxtn(:,1) &
      + b%bDf%uhxtn(:,2) )/2d0/b%dt


  end subroutine calcDerv
!!---------------------------End calcDerv--------------------------!!



!!---------------------------getVertVel----------------------------!!
  subroutine getVertVel(z,h,eta,hx,u,ux,uxx,uxxx,uhx,uhxx,uhxxx,&
    uxt,uhxt,uc,wc,pc)
  implicit none

    real(kind=8),intent(in)::z,h,eta,hx,u,ux,uxx,uxxx
    real(kind=8),intent(in)::uhx,uhxx,uhxxx,uxt,uhxt
    real(kind=8),intent(out)::uc,wc,pc

    uc = u - (0.5d0*h*uhxx - h*h*uxx/6d0) - (z*uhxx + 0.5d0*z*z*uxx)
    wc = -uhx - z*ux + z/2d0*( hx*uhxx + h*uhxxx ) &
      - z/6d0*( 2d0*h*hx*uxx + h*h*uxxx ) &
      + z*z/2d0*uhxxx + z*z*z/6d0*uxxx
    pc = rhoW * ( grav*( -z + eta ) + z*(uhxt + 0.5d0*z*uxt) )

  end subroutine getVertVel
!!-------------------------End getVertVel--------------------------!!
