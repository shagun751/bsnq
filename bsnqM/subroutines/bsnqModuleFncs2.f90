!!-----------------------------setMFree----------------------------!!
  subroutine setMFree(b)
  implicit none

    class(bsnqCase),intent(inout)::b
    integer(kind=C_K1)::nn,err,i,i2,j,k1
    integer(kind=C_K1),allocatable::neid(:),newrk(:)
    real(kind=C_K2)::cx,cy,rad,tmpr3,tmpr5
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
        write(9,'(" [ERR] No MFree at node ", I10)')i
        write(9,'(" [---] Cx, Cy ",2F15.6)')cx,cy
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

    !call b%bDf%init( b%npt )
    allocate(b%bDf( b%npt ))
    do i=1,b%npt
      b%bDf(i)%ux=0d0
      b%bDf(i)%uhx=0d0    
      b%bDf(i)%ux_tn=0d0
      b%bDf(i)%uhx_tn=0d0    
    enddo

    call system_clock(b%sysC(6))
    write(9,*)"[MSG] Done setMFree"
    write(9,'(" [TIM] ",F15.4)')1d0*(b%sysC(6)-b%sysC(5))/b%sysRate
    write(9,*)    

  end subroutine setMFree
!!---------------------------End setMFree--------------------------!!



!!-------------------------calcVertVelDerv-------------------------!!
  subroutine calcVertVelDerv(b)
  implicit none

    class(bsnqCase),intent(inout)::b
    integer(kind=C_K1)::i,i2,j,k,k2
    real(kind=C_K2)::tmpr1

    b%ur = b%tOb(0)%p / b%tOb(0)%tD
    b%vr = b%tOb(0)%q / b%tOb(0)%tD
    b%uhr = b%ur * b%dep
    b%vhr = b%vr * b%dep

    !$OMP PARALLEL DEFAULT(shared) &
    !$OMP   PRIVATE(i)
    !$OMP DO SCHEDULE(dynamic,100)    
    do i=1,b%npt
      b%bDf(i)%u = b%ur(i)
      ! Storing ! d(U)/dx and d(Uh)/dx at t(n-1), t(n-2)    
      b%bDf(i)%ux_tn(2) = b%bDf(i)%ux_tn(1)  
      b%bDf(i)%ux_tn(1) = b%bDf(i)%ux
      b%bDf(i)%uhx_tn(2) = b%bDf(i)%uhx_tn(1)    
      b%bDf(i)%uhx_tn(1) = b%bDf(i)%uhx
    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL
        
    !b%ur = dsin(b%cor(:,1))

    !! First derivative
    !$OMP PARALLEL DEFAULT(shared) &
    !$OMP   PRIVATE(i,i2,j,k,k2)
    !$OMP DO SCHEDULE(dynamic,100)    
    do i=1,b%npt
      i2 = b%pObf(i)%bsnqId
      k2 = b%pObf(i)%nn

      call calcGrad( k2, b%pObf(i)%phiDx, &
        b%ur( b%pObf(i)%neid ), b%bDf(i2)%ux, j )

      call calcGrad( k2, b%pObf(i)%phiDx, &
        b%uhr( b%pObf(i)%neid ), b%bDf(i2)%uhx, j )

      call calcGrad( k2, b%pObf(i)%phiDx, &
        b%dep( b%pObf(i)%neid ), b%bDf(i2)%hx, j )

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
        b%bDf(b%pObf(i)%neid)%ux, b%bDf(i2)%uxx, j )

      call calcGrad( k2, b%pObf(i)%phiDx, &
        b%bDf(b%pObf(i)%neid)%uhx, b%bDf(i2)%uhxx, j )

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
        b%bDf(b%pObf(i)%neid)%uxx, b%bDf(i2)%uxxx, j )

      call calcGrad( k2, b%pObf(i)%phiDx, &
        b%bDf(b%pObf(i)%neid)%uhxx, b%bDf(i2)%uhxxx, j )

    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL    

    tmpr1=2d0*b%dt
    !$OMP PARALLEL DEFAULT(shared) &
    !$OMP   PRIVATE(i)
    !$OMP DO SCHEDULE(dynamic,100)    
    do i=1,b%npt
      b%bDf(i)%uxt = ( 3d0*b%bDf(i)%ux - 4d0*b%bDf(i)%ux_tn(1) &
        + b%bDf(i)%ux_tn(2) )/tmpr1
      b%bDf(i)%uhxt = ( 3d0*b%bDf(i)%uhx - 4d0*b%bDf(i)%uhx_tn(1) &
        + b%bDf(i)%uhx_tn(2) )/tmpr1
    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL


  end subroutine calcVertVelDerv
!!-----------------------End calcVertVelDerv-----------------------!!



!!-------------------------findEleForLocXY-------------------------!!
  subroutine findEleForLocXY1(b,xin,yin,eleOut,ep,et)
  implicit none

    class(bsnqCase),intent(in)::b
    real(kind=C_K2),intent(in)::xin,yin
    integer(kind=C_K1),intent(out)::eleOut
    real(kind=C_K2),intent(out)::ep,et
    integer(kind=C_K1)::ind(3,2),p1,p2,dir,liel,li,nl(3)
    real(kind=C_K2)::xy(3,2),vec1(2),vec2(2),res

    ind(1,:)=(/ 1,2 /)
    ind(2,:)=(/ 2,3 /)
    ind(3,:)=(/ 3,1 /)

    do liel = 1, b%nele
      nl=b%conn(liel,1:3)
      do li=1,3
        xy(li,:)=b%cor(nl(li),:)
      enddo

      dir=1
      do li=1,3
        p1=ind(li,1)
        p2=ind(li,2)
        vec1(1)=xy(p2,1)-xy(p1,1)
        vec1(2)=xy(p2,2)-xy(p1,2)
        vec2(1)=xin-xy(p1,1)
        vec2(2)=yin-xy(p1,2)
        call vecCross2D(vec1(1),vec1(2),vec2(1),vec2(2),res)
        if(res.lt.0d0) dir=-1        
      enddo

      if(dir.eq.1)then
        eleOut=liel
        ep = ( b%invJ(liel,1)*(xin-xy(1,1)) ) &
          + ( b%invJ(liel,3)*(yin-xy(1,2)) )
        et = ( b%invJ(liel,2)*(xin-xy(1,1)) ) &
          + ( b%invJ(liel,4)*(yin-xy(1,2)) )
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
    integer(kind=C_K1)::ind(3,2),p1,p2,dir,pp,liel,li,nl(3)
    real(kind=C_K2)::xy(3,2),vec1(2),vec2(2),res,xp,yp

    ind(1,:)=(/ 1,2 /)
    ind(2,:)=(/ 2,3 /)
    ind(3,:)=(/ 3,1 /)


    !$OMP PARALLEL DEFAULT(shared) &
    !$OMP   PRIVATE(pp,xp,yp,liel,nl,xy,dir,li,p1,p2,vec1,vec2,res)
    !$OMP DO SCHEDULE(dynamic,10)
    do pp=1,np
      xp=xin(pp)
      yp=yin(pp)

      do liel = 1, b%nele
        nl=b%conn(liel,1:3)
        do li=1,3
          xy(li,:)=b%cor(nl(li),:)
        enddo

        dir=1
        do li=1,3
          p1=ind(li,1)
          p2=ind(li,2)
          vec1(1)=xy(p2,1)-xy(p1,1)
          vec1(2)=xy(p2,2)-xy(p1,2)
          vec2(1)=xp-xy(p1,1)
          vec2(2)=yp-xy(p1,2)
          call vecCross2D(vec1(1),vec1(2),vec2(1),vec2(2),res)
          if(res.lt.0d0) dir=-1        
        enddo

        if(dir.eq.1)then
          eleOut(pp)=liel
          ep(pp) = ( b%invJ(liel,1)*(xp-xy(1,1)) ) &
            + ( b%invJ(liel,3)*(yp-xy(1,2)) )
          et(pp) = ( b%invJ(liel,2)*(xp-xy(1,1)) ) &
            + ( b%invJ(liel,4)*(yp-xy(1,2)) )
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
!!-----------------------End findEleForLocXY-----------------------!!



!!----------------------------getVertVel---------------------------!!
  subroutine getVertVel(b,np,xin,yin,zin,uOut,vOut,wOut,pOut,&
    wrki,wrkr,err)
  implicit none

    class(bsnqCase),intent(in)::b
    integer(kind=C_K1),intent(in)::np
    real(kind=C_K2),intent(in)::xin(np),yin(np),zin(np)

    integer(kind=C_K1),intent(out)::wrki(np),err(np)
    real(kind=C_K2),intent(out)::uOut(np),vOut(np),wOut(np)
    real(kind=C_K2),intent(out)::pOut(np),wrkr(np,2)    

    integer(kind=C_K1)::nq(6),i,k
    real(kind=C_K2)::wei(6),hLoc,zLoc,etaLoc
    real(kind=C_K2)::uLoc,vLoc,wLoc,pLoc    
    type(vertVelDerv)::tmp

    !Note: pLoc here is pressure. It is not depth-integrated vel-x

    call b%findEleForLocXY2(np,xin,yin,wrki,wrkr(:,1),wrkr(:,2))

    !$OMP PARALLEL DEFAULT(shared) &
    !$OMP   PRIVATE(i, nq, wei, tmp, k, hLoc, etaLoc, zLoc, &
    !$OMP     uLoc, vLoc, wLoc, pLoc)
    !$OMP DO SCHEDULE(dynamic,10)
    do i = 1,np

      if(wrki(i).eq.-1)then
        err(i)=1
        uOut(i)=0d0
        vOut(i)=0d0
        wOut(i)=0d0
        pOut(i)=0d0     
        cycle   
      endif

      nq = b%conn(wrki(i),:)
      call fem_N6i(wrkr(i,1),wrkr(i,2),wei)      

      call tmp%initByInterp( 6, b%bDf(nq), wei, k )

      hLoc=0d0
      etaLoc=0d0      
      do k=1,6
        hLoc = hLoc + wei(k)*b%dep(nq(k))
        etaLoc = etaLoc + wei(k)*b%tOb(0)%tD(nq(k))
      enddo
      etaLoc = etaLoc - hLoc
      zLoc = zin(i) - hLoc

      call vertVelExp(zLoc, hLoc, etaLoc, tmp%hx,&
        tmp%u, tmp%ux, tmp%uxx, tmp%uxxx, &
        tmp%uhx, tmp%uhxx, tmp%uhxxx, &
        tmp%uxt, tmp%uhxt, uLoc, wLoc, pLoc)  

      vLoc=0d0

      err(i)=0
      uOut(i)=uLoc
      vOut(i)=vLoc
      wOut(i)=wLoc
      pOut(i)=pLoc

    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL    
    
  end subroutine getVertVel
!!--------------------------End getVertVel-------------------------!!



!!--------------------------testGetVertVel-------------------------!!
  subroutine testGetVertVel(b)
  implicit none

    class(bsnqCase),intent(in)::b
    
    integer,parameter::np=4
    integer(kind=C_K1)::wrki(np),err(np),i
    real(kind=C_K2)::xin(np),yin(np),zin(np),wrkr(np,2)
    real(kind=C_K2)::uOut(np),vOut(np),wOut(np),pOut(np)    
    
    xin(1)=15d0
    yin(1)=2.00d0
    zin(1)=0.20d0

    xin(2)=15d0
    yin(2)=2.00d0
    zin(2)=0.35d0

    xin(3)=16.02d0
    yin(3)=2.02d0
    zin(3)=0.20d0

    xin(4)=16.02d0
    yin(4)=2.02d0
    zin(4)=0.35d0

    call b%getVertVel(np,xin,yin,zin,uOut,vOut,wOut,pOut,&
      wrki,wrkr,err)

    write(120,'(F15.6)',advance='no')b%tOb(0)%rtm
    do i=1,np
      write(120,'(3F15.6)',advance='no')pOut(i),uOut(i),wOut(i)
    enddo
    write(120,*)

  end subroutine testGetVertVel
!!------------------------End testGetVertVel-----------------------!!



!!----------------------------locWvAng-----------------------------!!
  ! subroutine locWvAng(b)
  ! implicit none

  !   class(bsnqCase),intent(inout)::b

  !   integer(kind=C_K1)::i
  !   real(kind=C_K2)::p,q,pMag2
  !   real(kind=C_K2),parameter::velLowLimit2=1d-20

  !   !$OMP PARALLEL DEFAULT(shared) PRIVATE(i, p, q, pMag2)
  !   !$OMP DO SCHEDULE(dynamic,100)
  !   do i = 1, b%npt
  !     p = b%tOb(0)%p(i)
  !     q = b%tOb(0)%q(i)
  !     pMag2 = p**2 + q**2
  !     if(pMag2.lt.velLowLimit2) then 
  !       b%wvAng(i) = 0d0
  !     else  
  !       ! RESULT = ATAN2(Y, X) -pi to pi
  !       b%wvAng(i) = atan2(q, p)
  !     endif
  !   enddo
  !   !$OMP END DO NOWAIT
  !   !$OMP END PARALLEL    

  ! end subroutine locWvAng

  
  subroutine locWvAng(b, npt, eta6)
  implicit none

    class(bsnqCase),intent(inout)::b
    integer(kind=C_K1),intent(in)::npt
    real(kind=C_K2),intent(in)::eta6(npt)

    integer(kind=C_K1)::i, i2, nn, j
    real(kind=C_K2)::etx,ety

    !$OMP PARALLEL DEFAULT(shared) &
    !$OMP   PRIVATE(i, i2, nn, j, etx, ety)
    !$OMP DO SCHEDULE(dynamic,100)    
    do i = 1, b%npt
      i2 = b%pObf(i)%bsnqId
      nn = b%pObf(i)%nn

      call calcGrad( nn, b%pObf(i)%phiDx, &
        eta6( b%pObf(i)%neid ), etx, j )      

      call calcGrad( nn, b%pObf(i)%phiDy, &
        eta6( b%pObf(i)%neid ), ety, j )      
      
      ! RESULT = ATAN2(Y, X) -pi to pi
      b%wvAng(i2) = atan2(ety, etx)
    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

  end subroutine locWvAng
!!--------------------------End locWvAng---------------------------!!