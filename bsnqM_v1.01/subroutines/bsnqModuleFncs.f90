!!---------------------------dynaMatrices--------------------------!!
  subroutine dynaMatrices(b,tDr,ur,vr)
  implicit none

    class(bsnqCase),intent(inout)::b    
    real(kind=C_K2),intent(in)::tDr(b%npt),ur(b%npt),vr(b%npt)

    call matrixSet2(b%npl,b%npt,b%nele,b%conn,b%Sz,&
      b%ivl,b%ivq,b%linkl,b%linkq,b%invJ,b%dep,b%por,tDr,&
      ur,vr,b%gGx,b%gGy,b%gNAdv)

    !write(9,*)"[MSG] Done dynaMatrices"
    !write(9,*)

  end subroutine dynaMatrices
!!-------------------------End dynaMatrices------------------------!!



!!---------------------------statMatrices--------------------------!!
  subroutine statMatrices(b)
  implicit none

    class(bsnqCase),intent(inout)::b    

    call system_clock(b%sysC(5))
    call matrixSet1(b%npl,b%npt,b%nele,b%conn,b%Sz,b%ivl,b%ivq,&
      b%linkl,b%linkq,b%invJ,b%dep,b%por,b%mass1,b%mass2,&
      b%gBs1,b%gBs2,b%gBs3,b%gBs4,b%gCxF,b%gCyF,b%gDMat,&
      b%gBs5,b%gBs6)
    write(9,*)"[MSG] Done matrixSet1"

    call bndIntegral1(b%npl,b%npt,b%nele,b%nbnd,b%conn,b%mabnd,&
      b%Sz,b%ivl,b%ivq,b%linkl,b%linkq,b%invJ,b%bndS,b%dep,&
      b%gFBs1,b%gFBs2,b%gFBs3,b%gFBs4,b%gFW)
    write(9,*)"[MSG] Done bndIntegral1"

    b%massW=b%mass2
    b%massE=b%mass2
    b%gBs1=b%gBs1+b%gFBs1+b%mass1
    b%gBs2=b%gBs2+b%gFBs2
    b%gBs3=b%gBs3+b%gFBs3
    b%gBs4=b%gBs4+b%gFBs4+b%mass1
    b%gDMat=b%gDMat+b%gFW

    call diriBCMass(b%npl,b%npt,b%nbndp,b%bndP,b%bndPT,&
      b%Sz,b%ivl,b%ivq,b%linkl,b%linkq,b%bndPN,b%gBs1,b%gBs2,&
      b%gBs3,b%gBs4,b%massE)    
    write(9,*)"[MSG] Done dirichletBC"

    call b%CSRMatrices

    call b%destructR1

    call system_clock(b%sysC(6))
    write(9,*)"[MSG] Done statMatrices"
    write(9,'(" [TIM] ",F15.4)')1d0*(b%sysC(6)-b%sysC(5))/b%sysRate
    write(9,*)

  end subroutine statMatrices
!!-------------------------End statMatrices------------------------!!



!!----------------------------destructR1---------------------------!!
  subroutine destructR1(b)
  implicit none

    class(bsnqCase),intent(inout)::b    

    deallocate(b%gFBs1,b%gFBs2,b%gFBs3,b%gFBs4,b%gFW)
    deallocate(b%p2e,b%p2p)
    deallocate(b%aFull,b%ivf,b%linkf)
    deallocate(b%massW,b%massE)

  end subroutine destructR1
!!--------------------------End destructR1-------------------------!!



!!-----------------------------solveAll----------------------------!!
  subroutine solveAll(b,rkTime,pr,qr,pbpr,qbpr,er,&
    gXW,gXE,gXPQ,gRE,gRPQ)
  implicit none

    class(bsnqCase),intent(in)::b    
    real(kind=C_K2),intent(in)::pr(b%npt),qr(b%npt)
    real(kind=C_K2),intent(in)::pbpr(b%npt),qbpr(b%npt)
    real(kind=C_K2),intent(in)::er(b%npl),rkTime
    real(kind=C_DOUBLE),intent(out)::gXW(b%npl),gXE(b%npl),gRE(b%npl)
    real(kind=C_DOUBLE),intent(out)::gXPQ(2*b%npt),gRPQ(2*b%npt)
    real(kind=C_K2)::dt,absC

    !!  [Note] : 
    !!  Remeber everything is passed by reference
    !!  The vars pt1, er etc will be same as bq%pt1 and bq%et1
    !!  if you have passed them as the calling argument
    !!  The vars pt1, er etc will be same as bq%pt0 and bq%et0
    !!  if you have passed them as the calling argument
    !!  So modify the bq%et1 only after the entire computation
    !!  with its old values is done. Till then store it in bq%er

    dt=b%dt

    !!------------------solveW-----------------!!
    !gRE=0d0
    do i=1,b%npl
      k=(i-1)*b%ivl(0)
      tmpr1=0d0
      do j=1,b%ivl(i)
        k2=k+j
        i2=b%linkl(k2)
        tmpr1=tmpr1 + ( b%gDMat(k2)*er(i2) )
      enddo
      gRE(i)=tmpr1
    enddo

    gRE=gRE/b%rowMaxW

    call solveSys(b%npl,b%nnzl,b%ivsl,b%jvsl,b%gMW,gRE,gXW,&
      b%errLim,b%maxiter,i,tmpr1,j)
    write(9,301)'W',j,i,tmpr1
    !!----------------End solveW---------------!!    


    !!-----------------solveEta----------------!!
    !gRE=0d0
    do i=1,b%npl
      tmpr1=0d0
      tmpr2=0d0
      absC=-b%absC(i)

      ![3x6]
      k=(i-1)*b%ivq(0)      
      do j=1,b%ivq(i)
        k2=k+j
        i2=b%linkq(k2)
        tmpr1=tmpr1 + (b%gCxF(k2)*pr(i2)) &
            + (b%gCyF(k2)*qr(i2))
      enddo

      ![3x3]
      k=(i-1)*b%ivl(0)      
      do j=1,b%ivl(i)
        k2=k+j
        i2=b%linkl(k2)
        tmpr2=tmpr2 + ( absC*b%mass2(k2)*er(i2) )
      enddo

      gRE(i)=dt*( tmpr1 + tmpr2 )
    enddo

    call b%diriBCEtaDiff(gRE, b%tOb(0)%rtm, b%tOb(1)%rtm)
    
    gRE=gRE/b%rowMaxE

    !!  [Note] : 
    !!  Do not modify b%et0 and b%et1 yet, 
    !!  their old vals porbably were passed and calling arguments
    !!  and are required by the PQ equations
    call solveSys(b%npl,b%nnzl,b%ivsl,b%jvsl,b%gME,gRE,gXE,&
      b%errLim,b%maxiter,i,tmpr1,j)
    write(9,301)'Eta',j,i,tmpr1    
    call b%diriBCEtaDiff(gXE, b%tOb(0)%rtm, b%tOb(1)%rtm)
    !!---------------End solveEta--------------!!


    !!-----------------solvePQ-----------------!!    
    !gRPQ=0d0
    do i=1,b%npt      
      tmpr1=0d0
      tmpr2=0d0
      tmpr3=0d0
      tmpr4=0d0
      absC=-b%absC(i)
      
      ![6x3]
      k=(i-1)*b%ivl(0)
      do j=1,b%ivl(i) 
        k2=k+j
        i2=b%linkl(k2)        
        
        tmpr1=tmpr1 + (b%gGx(k2)*er(i2)) &
          + (b%gBs5(k2)*gXW(i2))

        tmpr2=tmpr2 + (b%gGy(k2)*er(i2)) &
          + (b%gBs6(k2)*gXW(i2))
      enddo

      ![6x6]
      k=(i-1)*b%ivq(0)
      do j=1,b%ivq(i)
        k2=k+j
        i2=b%linkq(k2)        
        
        tmpr3=tmpr3 + ( b%gNAdv(k2)*pbpr(i2) &
          + absC*b%mass1(k2)*pr(i2) )

        tmpr4=tmpr4 + ( b%gNAdv(k2)*qbpr(i2) &
          + absC*b%mass1(k2)*qr(i2) )
      enddo

      gRPQ(i)=dt*( tmpr1 + tmpr3 )
      gRPQ(b%npt+i)=dt*( tmpr2 + tmpr4 )
    enddo
    
    call b%diriBCPQDiff(gRPQ, b%tOb(0)%rtm, b%tOb(1)%rtm)

    gRPQ=gRPQ/b%rowMaxPQ
    
    call system_clock(sysC(1))
    call solveSys(2*b%npt,b%nnzf,b%ivsf,b%jvsf,b%gMPQ,gRPQ,&
      gXPQ,b%errLim,b%maxiter,i,tmpr1,j)
    write(9,301)'PQ',j,i,tmpr1
    call b%diriBCPQDiff(gXPQ, b%tOb(0)%rtm, b%tOb(1)%rtm)
    call system_clock(sysC(2))
    !!---------------End solvePQ---------------!!

    301 format('      |',a6,i10,i10,e15.4)

  end subroutine solveAll
!!---------------------------End solveAll--------------------------!!



!!-----------------------------initMat-----------------------------!!
  subroutine initMat(b)
  implicit none

    class(bsnqCase),intent(inout)::b

    i=b%npl
    j=b%npt
    i1=b%ivl(0)
    j1=b%ivq(0)

    allocate(b%gXE(i),b%gXPQ(2*j),b%gXW(i))
    allocate(b%por(j))
    allocate(b%ur(j),b%vr(j),b%pbpr(j),b%qbpr(j))

    allocate(b%massW(i1*i),b%massE(i1*i))
    allocate(b%mass1(j1*j),b%mass2(i1*i))    
    allocate(b%gBs1(j1*j),b%gBs2(j1*j))
    allocate(b%gBs3(j1*j),b%gBs4(j1*j))
    allocate(b%gCxF(j1*i),b%gCyF(j1*i),b%gDMat(i1*i))
    allocate(b%gBs5(i1*j),b%gBs6(i1*j),b%absC(j))
    allocate(b%gGx(i1*j),b%gGy(i1*j),b%gNAdv(j1*j))
    allocate(b%gFBs1(j1*j),b%gFBs2(j1*j))
    allocate(b%gFBs3(j1*j),b%gFBs4(j1*j),b%gFW(i1*i))
    allocate(b%aFull(b%ivf(0)*2*j))
    allocate(b%rowMaxW(i),b%rowMaxE(i),b%rowMaxPQ(2*j))
    allocate(b%gRE(i),b%gRPQ(2*j))
    allocate(b%etaMax(i),b%etaMin(i))

    b%Sz(1)=i1*i ![3x3] ![ivl(0) * npl]
    b%Sz(2)=j1*i ![3x6] ![ivq(0) * npl]
    b%Sz(3)=i1*j ![6x3] ![ivl(0) * npt]
    b%Sz(4)=j1*j ![6x6] ![ivq(0) * npt]
    b%por=1d0
    
    b%absC=0d0
    if(b%absOn)then
      do i=1,b%absOb(1)%N
        call b%absOb(i)%calcAbsC(b%npt,b%cor,b%absC)
      enddo  
    endif

    b%nTOb=2
    allocate(b%tOb(0:b%nTOb-1))
    do i=0,b%nTOb-1
      call b%tOb(i)%initBsnqVars(b%npl,b%npt)
    enddo

    b%nSOb=3
    allocate(b%sOb(b%nSOb))
    do i=1,b%nSOb
      call b%sOb(i)%initBsnqVars(b%npl,b%npt)
    enddo    

    b%tOb(0)%rtm=0d0
    b%tOb(0)%e=0d0
    b%tOb(0)%p=0d0
    b%tOb(0)%q=0d0
    b%tOb(0)%tD(1:b%npl)=b%dep(1:b%npl)+b%tOb(0)%e

    ! call solitIC(b%npl, b%npt, b%cor, &
    !   b%tOb(0)%e, b%tOb(0)%p, b%tOb(0)%q)
    
    b%etaMin=b%tOb(0)%e
    b%etaMax=b%tOb(0)%e
    call fillMidPoiVals(b%npl,b%npt,b%nele,b%conn,b%tOb(0)%tD)

    call paralution_init(b%nthrd)    

    call b%outputXML

    write(9,*)"[MSG] Done initMat"
    write(9,*)

  end subroutine initMat
!!---------------------------End initMat---------------------------!!



!!-----------------------------setRun------------------------------!!
  subroutine setRun(b)
  implicit none

    class(bsnqCase),intent(inout)::b

    !Input file open  
    bqtxt=trim(b%probname)//'.inp'
    inquire(file=trim(bqtxt),exist=ex)
    if(ex) then
      open(newunit=mf,file=trim(bqtxt))
    else
      write(9,*)"[ERR] Missing input file"
      stop
    endif

    write(9,'(" [MSG] setRun Unit = ",I10)')mf
    read(mf,*,end=81,err=81)bqtxt
    read(mf,*,end=81,err=81)b%resume
    read(mf,*,end=81,err=81)bqtxt
    read(mf,*,end=81,err=81)b%resumeFile  

    read(mf,*,end=81,err=81)bqtxt
    read(mf,*,end=81,err=81)b%dt     
    read(mf,*,end=81,err=81)bqtxt
    read(mf,*,end=81,err=81)b%endTime

    read(mf,*,end=81,err=81)bqtxt
    read(mf,*,end=81,err=81)b%errLim
    read(mf,*,end=81,err=81)bqtxt
    read(mf,*,end=81,err=81)b%maxIter

    read(mf,*,end=81,err=81)bqtxt
    read(mf,*,end=81,err=81)tmpr1
    b%fileOut=int(tmpr1/b%dt,4)
    read(mf,*,end=81,err=81)bqtxt
    read(mf,*,end=81,err=81)tmpr1
    b%resumeOut=int(tmpr1/b%dt,4)  

    read(mf,*,end=81,err=81)bqtxt
    read(mf,*,end=81,err=81)b%presOn      

    read(mf,*,end=81,err=81)bqtxt
    read(mf,*,end=81,err=81)b%nthrd
    read(mf,*,end=81,err=81)bqtxt
    read(mf,*,end=81,err=81)tmpi1
    if(tmpi1.gt.0) then
      allocate(b%probe(0:tmpi1))
      b%probe(0)=tmpi1
      read(mf,*,end=81,err=81)b%probe(1:b%probe(0))
    else
      allocate(b%probe(0:1))
      b%probe(0)=0
    endif

    read(mf,*,end=81,err=81)bqtxt
    read(mf,*,end=81,err=81)tmpr1,tmpr2,tmpr3
    read(mf,*,end=81,err=81)bqtxt
    read(mf,*,end=81,err=81)tmpr4,tmpr5,tmpr6
    ! waveType(T,d,H,X0,Y0,thDeg)
    b%wvIn=waveType(tmpr1,tmpr3,tmpr2,tmpr4,tmpr5,tmpr6)
    write(9,'(" [INF] ",3A15)')'T','L','d'
    write(9,'(" [---] ",3F15.6)')b%wvIn%T,b%wvIn%L,b%wvIn%d
    write(9,'(" [INF] ",A15)')'kh'
    write(9,'(" [---] ",F15.6)')b%wvIn%k*b%wvIn%d
    write(9,'(" [INF] At 2.25 WavePeriods")')
    write(9,'(" [---] ",3A15)')'Eta','P','Q'
    call b%wvIn%getEta(2.25d0*b%wvIn%T,b%wvIn%x0,b%wvIn%y0,tmpr1)
    call b%wvIn%getPQ(2.25d0*b%wvIn%T,b%wvIn%x0,b%wvIn%y0,tmpr2,tmpr3)
    write(9,'(" [---] ",3F15.6)')tmpr1,tmpr2,tmpr3

    read(mf,*,end=81,err=81)bqtxt
    read(mf,*,end=81,err=81)bqtxt
    read(mf,*,end=81,err=81)b%absOn
    read(mf,*,end=81,err=81)bqtxt
    read(mf,*,end=81,err=81)tmpi1    
    if(b%absOn.and.(tmpi1.le.0))then
      write(9,'(" [ERR] Improper number of absorbance layer")')
      stop
    endif
    if(tmpi1.gt.0)then
      allocate(b%absOb(tmpi1))
      do i=1,tmpi1
        read(mf,*,end=81,err=81)tmpi2,tmpr1,tmpr2,tmpr3
        b%absOb(i)=absTyp(tmpi1,tmpi2,tmpr1,tmpr2,tmpr3)
      enddo
    endif


    goto 82
    81 write(9,*) "[ERR] Check input file format"
    stop
    82 close(mf)

    write(9,*)"[MSG] Done setRun"
    write(9,*)
    
  end subroutine setRun
!!---------------------------End setRun----------------------------!!



!!----------------------------CSRMatrices--------------------------!!
  subroutine CSRMatrices(b)
  implicit none

    class(bsnqCase),intent(inout)::b    

    ! Full Matrice A
    do i=1,b%npt
      k=(i-1)*b%ivf(0)
      k2=(i-1)*b%ivq(0)
      b%aFull(k+1:k+b%ivq(i)-1)=b%gBs1(k2+1:k2+b%ivq(i)-1)
      b%aFull(k+b%ivq(i):k+b%ivf(i)-2)=b%gBs2(k2+1:k2+b%ivq(i)-1)
      b%aFull(k+b%ivf(i)-1)=b%gBs2(k2+b%ivq(i))
      b%aFull(k+b%ivf(i))=b%gBs1(k2+b%ivq(i))

      k=(i+b%npt-1)*b%ivf(0)
      k2=(i-1)*b%ivq(0)
      b%aFull(k+1:k+b%ivq(i)-1)=b%gBs3(k2+1:k2+b%ivq(i)-1)
      b%aFull(k+b%ivq(i):k+b%ivf(i)-2)=b%gBs4(k2+1:k2+b%ivq(i)-1)
      b%aFull(k+b%ivf(i)-1)=b%gBs3(k2+b%ivq(i))
      b%aFull(k+b%ivf(i))=b%gBs4(k2+b%ivq(i))
    enddo

    
    ! Normalising matrices
    b%rowMaxW=0d0
    b%rowMaxE=0d0
    b%rowMaxPQ=0d0
    do i=1,b%npl
      k=(i-1)*b%ivl(0)
      b%rowMaxW(i)=b%massW(k+b%ivl(i))
      b%rowMaxE(i)=b%massE(k+b%ivl(i))
      b%massW(k+1:k+b%ivl(i))=b%massW(k+1:k+b%ivl(i))/b%rowMaxW(i)
      b%massE(k+1:k+b%ivl(i))=b%massE(k+1:k+b%ivl(i))/b%rowMaxE(i)
      if((b%rowMaxW(i).eq.0d0).or.(b%rowMaxE(i).eq.0d0))then
        write(9,'(" [ERR] Check rowMaxE or rowMaxW, node",I10)')i
        write(9,'(" [---] rowMaxW",E15.6)')b%rowMaxW(i)
        write(9,'(" [---] rowMaxE",E15.6)')b%rowMaxE(i)
        stop
      endif
    enddo
    j1=2*b%npt
    j2=b%ivf(0)
    do i=1,j1
      k=(i-1)*j2
      k2=b%ivf(i)
      b%rowMaxPQ(i)=b%aFull(k+k2)
      b%aFull(k+1:k+k2)=b%aFull(k+1:k+k2)/b%rowMaxPQ(i)
      if(b%rowMaxPQ(i).eq.0d0)then
        write(9,'(" [ERR] Check rowMaxPQ, node, npt",2I10)')i,b%npt
        write(9,'(" [---] rowMaxPQ",E15.6)')b%rowMaxPQ(i)
        stop
      endif
    enddo

    
    ! Paralution CSR linear nodes
    b%nnzl=0
    do i=1,b%npl
      b%nnzl=b%nnzl+b%ivl(i)
    enddo

    allocate(b%ivsl(b%npl+1),b%jvsl(b%nnzl))
    allocate(b%gMW(b%nnzl),b%gME(b%nnzl))

    i2=0
    b%ivsl(1)=1
    do i=1,b%npl
      b%ivsl(i+1)=b%ivsl(i)+b%ivl(i)
      k=(i-1)*b%ivl(0)
      do j=1,b%ivl(i)
        k2=b%linkl(k+j)
        i2=i2+1
        b%jvsl(i2)=k2
        b%gMW(i2)=b%massW(k+j)
        b%gME(i2)=b%massE(k+j)
      enddo
    enddo
    if((i2.ne.b%nnzl).or.(b%ivsl(b%npl+1).ne.b%nnzl+1)) then
      write(9,*)'[ERR] CSR linear nnz not correct'
      stop
    endif


    ! Paralution CSR quadratic nodes
    b%nnzf=0
    do i=1,2*b%npt
      b%nnzf=b%nnzf+b%ivf(i)
    enddo

    allocate(b%ivsf(2*b%npt+1),b%jvsf(b%nnzf))
    allocate(b%gMPQ(b%nnzf))

    i2=0
    b%ivsf(1)=1
    do i=1,2*b%npt
      b%ivsf(i+1)=b%ivsf(i)+b%ivf(i)
      k=(i-1)*b%ivf(0)
      do j=1,b%ivf(i)
        k2=b%linkf(k+j)
        i2=i2+1
        b%jvsf(i2)=k2
        b%gMPQ(i2)=b%aFull(k+j)
      enddo
    enddo
    if((i2.ne.b%nnzf).or.(b%ivsf(2*b%npt+1).ne.b%nnzf+1)) then
      write(9,*)'[ERR] CSR quadratic nnz not correct'
      stop
    endif


    write(9,'(" [INF] Solve Lin ",2I10)')b%npl,b%nnzl
    write(9,'(" [INF] Solve Quad",2I10)')2*b%npt,b%nnzf
    write(9,*)"[MSG] Done CSRMatrices"

  end subroutine CSRMatrices
!!--------------------------End CSRMatrices------------------------!!



!!-----------------------------femInit-----------------------------!!
  subroutine femInit(b)
  implicit none

    class(bsnqCase),intent(inout)::b
    integer(kind=C_K1)::nbndpoi
    integer(kind=C_K1),allocatable::tempia(:,:)
    
    b%maxNePoi=maxNePoi
    b%maxNeEle=maxNeEle
    
    !!-----------------Generate Middle Points-----------------!!
    allocate(b%p2p(b%npt,b%maxNePoi),b%p2e(b%npt,b%maxNeEle))
    allocate(b%npoisur(b%npt,3))
    b%p2p=0
    b%npoisur=0

    call middleNode(b%npl,b%npq,b%npt,b%nele,b%nedg,&
      b%maxNePoi,b%cor,b%dep,b%conn,b%p2p,b%npoisur)
    if(b%npq.ne.b%nedg) then
      write(9,*) "[ERR] Initial npoint assumption insufficiant"
      stop
    endif
    b%p2p=0
    b%npoisur=0

    !Debug Comments
    ! do i=1,npoint
    !   write(9,*)i,":",(coord(i,j),j=1,2),depth(i)
    ! enddo
    ! do i=1,nelem
    !   write(9,*)i,":",(conn(i,j),j=1,6)
    ! enddo
    write(9,*)"[MSG] Middle points generation done"
    write(9,'(" [INF] ",4A12)')"LinNode","QuadNode","TotNode",&
      "nEdges"
    write(9,'(" [---] ",4I12)')b%npl,b%npq,b%npt,b%nedg

    !Boundary sides middle nodes update
    !middle nodes stored in mabnd(i,5)
    call bndMiddleNodes(b%npl,b%npt,b%nele,b%nbnd,&
      b%conn,b%mabnd)
    write(9,*) "[MSG] Boundary sides middle point done"  
    !Debug Comments
    ! do i=1,nbnd
    !   write(9,*)i,":",(mabnd(i,j),j=1,6)
    ! enddo

    !Storing boundary nodes in bndNode
    nbndpoi=0
    allocate(tempia(b%npt,2))
    tempia=0
    do i=1,b%nbnd
      nq(1)=b%mabnd(i,1)
      nq(2)=b%mabnd(i,2)
      nq(3)=b%mabnd(i,5)
      tmpi4=b%mabnd(i,4)
      do j=1,3
        do k=1,nbndpoi
          if(tempia(k,1).eq.nq(j)) goto 31
        enddo
        nbndpoi=nbndpoi+1
        tempia(nbndpoi,1)=nq(j)
        tempia(nbndpoi,2)=tmpi4
        31 continue
        !Preferance for bndtype 14 over 12 or 13
        if((tempia(k,2).ne.14).and.(tmpi4.eq.14)) then
          tempia(k,2)=14
        end if
        !Preferance for bndtype 11 - higher
        if((tempia(k,2).ne.11).and.(tmpi4.eq.11)) then
          tempia(k,2)=11
        end if
        ! if((tempia(k,2).eq.12).and.(tmpi4.eq.13)) then
        !   tempia(k,2)=13
        ! end if
      enddo
    enddo
    allocate(b%bndP(nbndpoi),b%bndPT(nbndpoi))
    b%bndP=tempia(1:nbndpoi,1)
    b%bndPT=tempia(1:nbndpoi,2)
    b%nbndp=nbndpoi
    deallocate(tempia)
    allocate(tempia(nbndpoi,2))
    tempia(:,1)=b%bndP
    tempia(:,2)=b%bndPT
    b%bndPT=0
    call mergeSort(nbndpoi,nbndpoi,b%bndP)
    do i=1,nbndpoi
      tmpi1=b%bndP(i)
      do j=1,nbndpoi
        if(tmpi1.eq.tempia(j,1)) goto 35
      enddo
      35 b%bndPT(i)=tempia(j,2)
    enddo
    deallocate(tempia)
    ! write(9,*)nbndpoi
    ! do i=1,nbndpoi
    !   k=bndNode(i)
    !   write(9,*)i,":",bndNodeType(i),coord(k,1),coord(k,2)
    ! enddo
  !!-----------------End Generate Middle Points---------------!!

  !!----------------------Node Connectivity-------------------!!
    call nodeConnVSR(b%npt,b%nele,b%maxNePoi,b%maxNeEle,&
      b%conn,b%p2p,b%p2e,b%npoisur)
    do i=1,b%npt
      call mergeSort(b%maxNePoi,b%npoisur(i,1),b%p2p(i,:))
    enddo

    !Finding number of linear element nbhs
    do i=1,b%npt
      do j=1,b%npoisur(i,1)
        if(b%p2p(i,j).gt.b%npl) exit
        b%npoisur(i,3)=b%npoisur(i,3)+1
      enddo
    enddo

    !Debug Comments
    ! do i=1,npoint
    !   write(9,*) i,":",(poi2poi(i,j),j=1,maxNePoi)
    ! enddo
    ! do i=1,npoint
    !   write(9,*) i,":",(poi2ele(i,j),j=1,maxNeEle)
    ! enddo
    ! do i=1,npoint
    !   write(9,*) i,":",(npoisur(i,j),j=1,3)
    ! enddo
    write(9,*)"[MSG] Node Connectivity Done"
  !!--------------------End Node Connectivity-----------------!!

  !!--------------------------VSR Matrices--------------------!!
    !Storage allocations
    allocate(b%ivl(0:b%npt),b%ivq(0:b%npt))  
    b%ivl(0)=maxval(b%npoisur(:,3))+1
    b%ivq(0)=maxval(b%npoisur(:,1))+1
    allocate(b%linkl(b%ivl(0)*b%npt),b%linkq(b%ivq(0)*b%npt))
    b%linkl=0
    b%linkq=0

    !IV matrix linear and quadratic
    do i=1,b%npl
      b%ivl(i)=b%npoisur(i,3)+1
    enddo
    do i=b%npl+1,b%npt
      b%ivl(i)=b%npoisur(i,3)
    enddo
    ! write(9,*) ivl

    do i=1,b%npt
      b%ivq(i)=b%npoisur(i,1)+1
    enddo
    ! write(9,*) ivq

    !link matrix linear
    do i=1,b%npl
      k=(i-1)*b%ivl(0)
      do j=1,b%npoisur(i,3)
        b%linkl(k+j)=b%p2p(i,j)
      enddo
      b%linkl(k+b%ivl(i))=i
    enddo
    do i=b%npl+1,b%npt
      k=(i-1)*b%ivl(0)
      do j=1,b%npoisur(i,3)
        b%linkl(k+j)=b%p2p(i,j)
      enddo    
    enddo
    ! do i=1,npoint
    !   k=(i-1)*ivl(0)
    !   write(9,*) i,":",(linkl(k+j),j=1,ivl(i))
    ! enddo

    !link matrix quadratic
    do i=1,b%npt
      k=(i-1)*b%ivq(0)
      do j=1,b%npoisur(i,1)
        b%linkq(k+j)=b%p2p(i,j)
      enddo
      b%linkq(k+b%ivq(i))=i
    enddo
    ! write(9,*)linkq

    ! Full Matrices ivFull and linkf
    allocate(b%ivf(0:2*b%npt))
    b%ivf(0)=2*b%ivq(0)
    allocate(b%linkf(b%ivf(0)*2*b%npt))
    do i=1,b%npt
      b%ivf(i)=2*b%ivq(i)
      b%ivf(b%npt+i)=2*b%ivq(i)
    enddo
    
    do i=1,b%npt
      k=(i-1)*b%ivf(0)
      b%linkf(k+1:k+b%npoisur(i,1))=b%p2p(i,1:b%npoisur(i,1))
      b%linkf(k+b%npoisur(i,1)+1:k+(2*b%npoisur(i,1)))=b%npt &
        + b%p2p(i,1:b%npoisur(i,1))
      b%linkf(k+b%ivf(i)-1)=b%npt+i
      b%linkf(k+b%ivf(i))=i

      k=(i+b%npt-1)*b%ivf(0)
      b%linkf(k+1:k+b%npoisur(i,1))=b%p2p(i,1:b%npoisur(i,1))
      b%linkf(k+b%npoisur(i,1)+1:k+(2*b%npoisur(i,1)))=b%npt &
        + b%p2p(i,1:b%npoisur(i,1))
      b%linkf(k+b%ivf(i)-1)=i
      b%linkf(k+b%ivf(i))=b%npt+i
    enddo

    write(9,*)"[MSG] VSR storage matrices done"
  !!-----------------------End VSR Matrices-------------------!!

  !!---------------------Jacobian and Normals-----------------!!
    allocate(b%invJ(b%nele,5),b%bndS(b%nbnd,3))
    call jacbInvLin(b%npt,b%nele,b%conn,b%cor,b%invJ)      
    ! Priniting area using linear and quad jacb
    ! do iel=1,nelem
    !   tmpr1=invJ(iel,5)/2d0
    !   tmpr2=0d0
    !   do i=1,nGP
    !     tmpr2=tmpr2+gpW(i)*jacb(iel)%D(i)
    !   enddo
    !   !if(abs(tmpr2-tmpr1).lt.1e-10)cycle
    !   write(*,'(i10,3f15.8)')iel,tmpr1,tmpr2,abs(tmpr2-tmpr1)    
    ! enddo  
    ! stop

    ! Boundary side normals
    call bndSideInfo(b%npl,b%npt,b%nele,b%nbnd,b%cor,b%mabnd,&
      b%bndS)
    ! do i=1,nbnd
    !   write(9,*)i,":",bndSide(i,:),mabnd(i,1:5)
    ! enddo  

    !Boundary Node Normal
    allocate(b%bndPN(b%npt,2))
    call bndNodeNormal(b%npl,b%npt,b%nele,b%nbnd,b%cor,b%mabnd,&
      b%bndS,b%bndPN)
    ! do i=1,npoint
    !   write(9,*)i,":",bndNormal(i,1),bndNormal(i,2)
    ! enddo
    write(9,*)"[MSG] Done femInit"
    write(9,*)
  !!----------------End Jacobian and Normals------------------!!

  end subroutine femInit
!!---------------------------End femInit---------------------------!!



!!-----------------------------meshRead----------------------------!!
  subroutine meshRead(b)
  implicit none

    class(bsnqCase),intent(inout)::b

    call system_clock(b%sysC(9),b%sysC(10))
    b%sysRate=real(b%sysC(10),C_K2)
    write(9,'(" [MSG] sysRate = ",F20.4)')b%sysRate    

    !Opening mesh file  
    bqtxt=trim(b%probname)//'.plt'
    inquire(file=trim(bqtxt),exist=ex)
    if(ex) then
      open(newunit=mf,file=trim(bqtxt))    
    else
      write(*,*)"[ERR] Missing mesh file"
      stop
    endif

    !!-------------------------Mesh Type 0--------------------!!
    write(9,'(" [MSG] meshRead Unit = ",I10)')mf
    write(9,'(" [MSG] C_K1, C_K2 = ",2I10)')C_K1,C_K2
    read(mf,*,end=11,err=11)bqtxt
    read(mf,*,end=11,err=11)b%nele,b%npl,b%nedg
    read(mf,*,end=11,err=11)bqtxt
    read(mf,*,end=11,err=11)b%nbnd,b%nbndtyp
    write(9,'(3a15)')"Elements","Linear Nodes","Edges"
    write(9,'(3i15)')b%nele,b%npl,b%nedg
    write(9,'(2a15)')"Bnd","BndTypes"
    write(9,'(2i15)')b%nbnd,b%nbndtyp

    !Assumption regarding number of quad Nodes
    b%npq=b%nedg
    b%npt=b%npl+b%npq

    !Nodes Read
    allocate(b%cor(b%npt,2))
    b%cor=-999  
    read(mf,*,end=11,err=11)bqtxt
    do i=1,b%npl
      read(mf,*,end=11,err=11)b%cor(i,1),b%cor(i,2)
    enddo
    write(9,*)"[MSG] Nodes read done"
    !Debug Comments
    ! do i=1,npoinl
    !   write(9,*)i,":",(coord(i,j),j=1,2),depth(i)
    ! enddo

    !Elements Read
    allocate(b%conn(b%nele,6))
    b%conn=0
    read(mf,*,end=11,err=11)bqtxt
    do i=1,b%nele
      read(mf,*,end=11,err=11)b%conn(i,3),&
        b%conn(i,1),b%conn(i,2)    
    enddo
    write(9,*)"[MSG] Elements read done"
    !Debug Comments
    ! do i=1,nelem
    !   write(9,*)i,":",(conn(i,j),j=1,6)
    ! enddo

    !Boundaries Read
    allocate(b%mabnd(b%nbnd,6))
    b%mabnd=0
    k=0;
    do j=1,b%nbndtyp
      read(mf,*,end=11,err=11)bqtxt
      read(mf,*,end=11,err=11)tmpi1,tmpi2  
      do i=k+1,k+tmpi2
        read(mf,*,end=11,err=11)b%mabnd(i,1:3)
        b%mabnd(i,4)=tmpi1      
      enddo
      k=k+tmpi2    
    enddo
    if(k.ne.b%nbnd) goto 14
    write(9,*) "[MSG] Boundaries read done" 
    ! ! Debug Comments 
    ! do i=1,nbnd
    !   write(9,*)i,":",(mabnd(i,j),j=1,6)
    ! enddo

    !Depth Read
    allocate(b%dep(b%npt))
    b%dep=-999
    read(mf,*,end=11,err=11)bqtxt
    do i=1,b%npl
      read(mf,*,end=11,err=11)b%dep(i)
    enddo
    write(9,*) "[MSG] Depth read done"
    !Debug Comments
    ! do i=1,npoinl
    !   write(9,*)i,depth(i)
    ! enddo


    goto 12
    !!-----------------------End Mesh Type 0------------------!!

    11 write(9,*) "[ERR] Check mesh file format"
    stop
    13 write(9,*) "[ERR] hex2dec error"
    stop
    14 write(9,*) "[ERR] Number of boundaries mismatch"
    stop
    12 close(mf)
    write(9,*)"[MSG] Done meshRead "
    write(9,*)
  end subroutine meshRead
!!---------------------------End meshRead--------------------------!!