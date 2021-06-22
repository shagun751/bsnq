!!-------------------------airyWaveModule--------------------------!!
module airyWaveModule
use bsnqGlobVars
implicit none

  type, public :: airyType    
    real(kind=C_K2)::T,d,H,L,w
    real(kind=C_K2)::x0,y0
    real(kind=C_K2)::k,kx,ky
    real(kind=C_K2)::thDeg,thRad,csth,snth      
    !To make cos signal start from eta=0 with crest
    real(kind=C_K2)::phi0 
  contains
    procedure ::  getEta
    procedure ::  getPQ
  end type airyType

  interface airyType
     procedure :: waveLenCalc
  end interface airyType

contains

!!-------------------------waveLenCalc-------------------------!!
  type(airyType) function waveLenCalc(inT,inD,inH,inX0,inY0,inThDeg)
  implicit none

    !! WaveLen using dispersion relation from Airy wave-theory

    integer(kind=C_K1)::iterMax,i
    real(kind=C_K2)::l0,newl,oldl,x,errLim
    real(kind=C_K2),intent(in)::inT,inD,inH,inX0,inY0,inThDeg


    waveLenCalc%T=inT
    waveLenCalc%d=inD
    waveLenCalc%H=inH
    waveLenCalc%x0=inX0
    waveLenCalc%y0=inY0
    waveLenCalc%thDeg=inThDeg
    waveLenCalc%w=2d0*pi/waveLenCalc%T
    waveLenCalc%thRad=waveLenCalc%thDeg*pi/180d0
    waveLenCalc%csth=dcos(waveLenCalc%thRad)
    waveLenCalc%snth=dsin(waveLenCalc%thRad)

    !To make cos signal start from eta=0 with crest
    waveLenCalc%phi0 = pi/2d0

    iterMax=50000
    errLim=1d-6
    l0 = (grav/2d0/pi)*(waveLenCalc%T)**2
    oldl = l0  
    do i = 1,iterMax
      newl = l0*dtanh(2d0*pi*(waveLenCalc%d)/oldl)
      !if(mod(i,100).eq.0) write(*,*)i,oldl,newl    
      x = abs(newl-oldl)
      if (x.le.errLim) then
        waveLenCalc%L = newl
        exit
      else
        oldl = newl
      end if
    end do

    if(i.ge.iterMax) then
      write(*,*)
      write(*,*)"[ERR] waveCalculator Error waveL",waveLenCalc%L
      write(*,*)
      waveLenCalc%L=-999
      waveLenCalc%k=0d0
      return
    endif

    waveLenCalc%k=2d0*pi/waveLenCalc%L
    waveLenCalc%kx=waveLenCalc%k*waveLenCalc%csth
    waveLenCalc%ky=waveLenCalc%k*waveLenCalc%snth

    write(9,'(" [INF] Airy. Vel Wheeler Stretching")')
    write(9,'(" [---] Vel Integrated [-h 0] and multiplied by (h+eta)/h ")')
    write(9,'(" [---] because of Wheeler stretching. Ref doc log_bsnqM_vAlgo.md")')
    write(9,*)

  end function waveLenCalc
!!-----------------------End waveLenCalc-----------------------!!



!!---------------------------getEta----------------------------!!
  subroutine getEta(b,rTime,x,y,eta)
  implicit none

    class(airyType),intent(in)::b

    !integer(kind=C_K1)::

    real(kind=C_K2),intent(in)::rTime,x,y
    real(kind=C_K2),intent(out)::eta
    real(kind=C_K2)::dx,dy

    dx=(x-b%x0)
    dy=(y-b%y0)
    eta=b%H/2d0 * dcos(b%kx*dx + b%ky*dy - b%w*rTime + b%phi0)

  end subroutine getEta
!!-------------------------End getEta--------------------------!!



!!----------------------------getPQ----------------------------!!
  subroutine getPQ(b,rTime,x,y,eta,p,q)
  implicit none

    class(airyType),intent(in)::b

    !integer(kind=C_K1)::

    real(kind=C_K2),intent(in)::rTime,x,y,eta
    real(kind=C_K2),intent(out)::p,q
    real(kind=C_K2)::dx, dy, phi, pn 
    
    dx=(x-b%x0)
    dy=(y-b%y0)
    phi = b%kx*dx + b%ky*dy - b%w*rTime + b%phi0

    ! Integrating from -h to 0
    pn = b%H/2d0 * b%w / b%k * dcos(phi) 

    ! Using Wheeler stretching
    pn = pn * (b%d + eta)/b%d

    p = pn * b%csth 
    q = pn * b%snth 

  end subroutine getPQ
!!--------------------------End getPQ--------------------------!!


end module airyWaveModule
!!-----------------------End airyWaveModule------------------------!!







!!------------------------stokes2WaveModule------------------------!!
module stokes2WaveModule
use bsnqGlobVars
implicit none

  type, public :: stokes2Type    
    real(kind=C_K2)::T,d,H,L,w
    real(kind=C_K2)::x0,y0
    real(kind=C_K2)::k,kx,ky
    real(kind=C_K2)::thDeg,thRad,csth,snth      
    !To make cos signal start from eta=0 with crest
    real(kind=C_K2)::phi0 
  contains
    procedure ::  getEta
    procedure ::  getPQ
  end type stokes2Type

  interface stokes2Type
     procedure :: waveLenCalc
  end interface stokes2Type

contains

!!-------------------------waveLenCalc-------------------------!!
  type(stokes2Type) function waveLenCalc(inT,inD,inH,inX0,inY0,inThDeg)
  implicit none

    !! WaveLen using dispersion relation from Airy wave-theory

    integer(kind=C_K1)::iterMax,i
    real(kind=C_K2)::l0,newl,oldl,x,errLim
    real(kind=C_K2)::qa, qb, phi0, kh
    real(kind=C_K2),intent(in)::inT,inD,inH,inX0,inY0,inThDeg


    waveLenCalc%T=inT
    waveLenCalc%d=inD
    waveLenCalc%H=inH
    waveLenCalc%x0=inX0
    waveLenCalc%y0=inY0
    waveLenCalc%thDeg=inThDeg
    waveLenCalc%w=2d0*pi/waveLenCalc%T
    waveLenCalc%thRad=waveLenCalc%thDeg*pi/180d0
    waveLenCalc%csth=dcos(waveLenCalc%thRad)
    waveLenCalc%snth=dsin(waveLenCalc%thRad)

    iterMax=50000
    errLim=1d-6
    l0 = (grav/2d0/pi)*(waveLenCalc%T)**2
    oldl = l0  
    do i = 1,iterMax
      newl = l0*dtanh(2d0*pi*(waveLenCalc%d)/oldl)
      !if(mod(i,100).eq.0) write(*,*)i,oldl,newl    
      x = abs(newl-oldl)
      if (x.le.errLim) then
        waveLenCalc%L = newl
        exit
      else
        oldl = newl
      end if
    end do

    if(i.ge.iterMax) then
      write(*,*)
      write(*,*)"[ERR] waveCalculator Error waveL",waveLenCalc%L
      write(*,*)
      waveLenCalc%L=-999
      waveLenCalc%k=0d0
      return
    endif

    waveLenCalc%k=2d0*pi/waveLenCalc%L
    waveLenCalc%kx=waveLenCalc%k*waveLenCalc%csth
    waveLenCalc%ky=waveLenCalc%k*waveLenCalc%snth

    ! Calculation phi0 to make the wave start from
    ! eta=0 followed by a crest
    kh = waveLenCalc%k * waveLenCalc%d
    qa = ( waveLenCalc%H/2d0  )
    qb = ( waveLenCalc%H**2 * waveLenCalc%k/16d0 &
      * dcosh(kh)/(dsinh(kh)**3) * (2d0 + dcosh(2*kh)) )

    phi0 = (-qa + dsqrt(qa*qa + 8d0*qb*qb))/(4d0*qb)
    waveLenCalc%phi0 = acos(phi0)
    write(9,'(" [INF] Stokes2. Vel Integrated [-h eta]")')
    write(9,'(" [INF] Stokes2 phase-0 radians", 2F15.6)') &
      waveLenCalc%phi0
    write(9,*)

  end function waveLenCalc
!!-----------------------End waveLenCalc-----------------------!!



!!---------------------------getEta----------------------------!!
  subroutine getEta(b,rTime,x,y,eta)
  implicit none

    class(stokes2Type),intent(in)::b

    !integer(kind=C_K1)::

    real(kind=C_K2),intent(in)::rTime,x,y
    real(kind=C_K2),intent(out)::eta
    real(kind=C_K2)::dx, dy, phi, kh

    dx=(x-b%x0)
    dy=(y-b%y0)
    phi = b%kx*dx + b%ky*dy - b%w*rTime + b%phi0
    kh = b%k * b%d
    eta = ( b%H/2d0 * dcos(phi) ) + &
      ( b%H**2 * b%k/16d0 * dcosh(kh)/(dsinh(kh)**3) &
        * (2d0 + dcosh(2*kh)) * dcos(2d0*phi) )

  end subroutine getEta
!!-------------------------End getEta--------------------------!!



!!----------------------------getPQ----------------------------!!
  subroutine getPQ(b,rTime,x,y,eta,p,q)
  implicit none

    class(stokes2Type),intent(in)::b

    !integer(kind=C_K1)::

    real(kind=C_K2),intent(in)::rTime,x,y,eta
    real(kind=C_K2),intent(out)::p, q
    real(kind=C_K2)::dx, dy, phi, kh, pn

    
    dx=(x-b%x0)
    dy=(y-b%y0)
    phi = b%kx*dx + b%ky*dy - b%w*rTime + b%phi0
    kh = b%k * b%d

    ! No wheeler stretching, Integrating -h to eta
    ! First order
    pn = grav * b%H/2d0 / b%w / dcosh(kh) &
      * dsinh(kh + b%k*eta) * dcos(phi)
    ! Second order
    pn = pn + ( 3d0/32d0 * b%H**2 * b%w / (dsinh(kh)**4) &
      * dsinh(2d0*b%k*(b%d+eta)) ) * dcos(2d0*phi)    

    ! Wheeler stretching and Integrating -h to 0
    ! ! First order
    ! pn = b%H/2d0 * b%w / b%k * dcos(phi)
    ! ! Second order
    ! pn = pn + ( 3d0/16d0 * b%H**2 * b%w / dtanh(kh)  &
    !   / (dsinh(kh)**2) ) * dcos(2d0*phi)    
    ! pn = pn * (b%d + eta)/b%d

    p = pn * b%csth

    q = pn * b%snth      

  end subroutine getPQ
!!--------------------------End getPQ--------------------------!!

end module stokes2WaveModule
!!----------------------End stokes2WaveModule----------------------!!







!!--------------------------outAbsModule---------------------------!!
module outAbsModule
use bsnqGlobVars
implicit none

  type, public :: absTyp
    integer(kind=C_K1)::N,typ
    real(kind=C_K2)::x0,l,c1,wT
  contains
    procedure ::  calcAbsC
  end type absTyp

  interface absTyp
    procedure ::  initAbsTyp    
  end interface absTyp

contains

!!--------------------------initAbsTyp-------------------------!!
  function initAbsTyp(N,typ,x0,l,wT) result(b)
  implicit none
    
    integer(kind=C_K1),intent(in)::N,typ
    real(kind=C_K2),intent(in)::x0,l,wT    
    type(absTyp)::b    

    b%N=N
    b%typ=typ
    b%x0=x0
    b%l=l
    b%wT=wT
    b%c1=30d0/wT/(dexp(1d0)-1d0)    

  end function initAbsTyp
!!------------------------End initAbsTyp-----------------------!!



!!---------------------------calcAbsC--------------------------!!
  subroutine calcAbsC(b,npt,cor,absC)
  implicit none

    class(absTyp),intent(in)::b
    integer(kind=C_K1),intent(in)::npt
    integer(kind=C_K1)::i,ix

    real(kind=C_K2),intent(in)::cor(npt,2)
    real(kind=C_K2),intent(inout)::absC(npt)
    real(kind=C_K2)::x,dc,dx

    select case(b%typ)
      case(1)
        ix=2
        dc=1d0

      case(2)
        ix=1
        dc=-1d0

      case(3)
        ix=2
        dc=-1d0

      case(4)
        ix=1
        dc=1d0

    end select

    do i=1,npt
      dx=dc*(cor(i,ix)-b%x0)/b%l
      if(dx.gt.0d0)then
        dx = min(dx, 1d0)
        absC(i)=b%c1*(dexp(dx**2)-1d0)
      endif
    enddo

  end subroutine calcAbsC
!!-------------------------End calcAbsC------------------------!!

end module outAbsModule
!!------------------------End outAbsModule-------------------------!!







!!-------------------------waveFileModule--------------------------!!
module waveFileModule
use bsnqGlobVars
implicit none

  type, public :: wvFileType        
    character(len=256)::fileName
    integer(kind=C_K1)::numP,posI
    real(kind=C_K2),allocatable::data(:,:), datadtt(:,:)
    real(kind=C_K2)::rampt0, rampt1, rampdt
  contains
    procedure ::  chkWaveFileInput
    procedure ::  initWaveFile
    procedure ::  initAiryFile
    procedure ::  initStokes2File
    procedure ::  getEta
    procedure ::  getPQ
  end type wvFileType  

contains

!!------------------------initWaveFile-------------------------!!
  subroutine initWaveFile(b, rampt0, rampt1)
  implicit none

    class(wvFileType),intent(inout)::b
    real(kind=C_K2),intent(in)::rampt0, rampt1
    integer(kind=C_K1)::mf,i
    real(kind=C_K2)::tmpra(4), ramptw, t
    logical::ex

    inquire(file=trim(b%fileName),exist=ex)
    if(ex) then
      open(newunit=mf,file=trim(b%fileName))
    else
      write(9,*)"[ERR] Missing wave input file"
      stop
    endif

    b%numP=0
    do while(.true.)
      read(mf,*,end=11,err=12)tmpra
      b%numP=b%numP+1
      cycle
      11 exit
      12 write(9,*)"[ERR] Check wave file format"
      stop
    enddo
    write(9,'(" [INF] Number of wave file point = ",i10)')b%numP
    close(mf)

    b%rampt0 = rampt0
    b%rampt1 = rampt1
    b%rampdt = (rampt1 - rampt0)

    allocate(b%data(b%numP,4), b%datadtt(b%numP,4))
    open(newunit=mf,file=trim(b%fileName))    
    do i=1,b%numP
      read(mf,*,end=21,err=21)b%data(i,1:4) 

      t = b%data(i,1)

      !Time Ramp
      ramptw = 1d0
      if(t.le.b%rampt0)then
        ramptw = 0d0
      else if(t.ge.b%rampt1)then
        ramptw = 1d0
      else if((t.gt.b%rampt0) .and. (t.lt.b%rampt1))then
        ramptw = 0.5d0*(1d0 - dcos(pi*(t - b%rampt0)/b%rampdt) )
      endif

      b%data(i,2:4) = b%data(i,2:4) * ramptw
    enddo
    b%datadtt(:,1) = b%data(:,1)
    write(9,*)"[INF] Done wave file read"

    call setCubicSpline(b%numP, 4, b%data, b%datadtt)
    b%posI=1

    return
    21 write(9,*)"[ERR] Check wave file format"
    stop
  end subroutine initWaveFile
!!----------------------End initWaveFile-----------------------!!



!!------------------------initAiryFile-------------------------!!
  subroutine initAiryFile(b,dt,inTotT,inT,inD,inH,inAng,&
    rampt0, rampt1)
  use airyWaveModule
  implicit none

    class(wvFileType),intent(inout)::b
    integer(kind=C_K1)::i
    real(kind=C_K2),intent(in)::inTotT,inT,inD,inH,dt,inAng
    real(kind=C_K2),intent(in)::rampt0, rampt1
    real(kind=C_K2)::t, eta, p, q, totT, ramptw
    type(airyType)::wv    

    ! airyType(T,d,H,X0,Y0,thDeg)
    wv=airyType(inT,inD,inH,0d0,0d0,inAng)

    write(9,'(" [INF] ",3A15)')'T','L','d'
    write(9,'(" [---] ",3F15.6)')wv%T,wv%L,wv%d
    write(9,'(" [INF] ",A15)')'kh'
    write(9,'(" [---] ",F15.6)')wv%k*wv%d
    write(9,'(" [INF] At 2.25 WavePeriods")')
    write(9,'(" [---] ",3A15)')'Eta','P','Q'
    call wv%getEta(2.25d0*wv%T,wv%x0,wv%y0,eta)
    call wv%getPQ(2.25d0*wv%T,wv%x0,wv%y0,eta,p,q)
    write(9,'(" [---] ",3F15.6)')eta,p,q

    totT=1.2d0*inTotT
    b%numP=floor(totT/dt)+2
    b%rampt0 = rampt0
    b%rampt1 = rampt1
    b%rampdt = (rampt1 - rampt0)

    allocate(b%data(b%numP,4), b%datadtt(b%numP,4))
    do i=0,b%numP-1
      t=i*dt

      !Time Ramp
      ramptw = 1d0
      if(t.le.b%rampt0)then
        ramptw = 0d0
      else if(t.ge.b%rampt1)then
        ramptw = 1d0
      else if((t.gt.b%rampt0) .and. (t.lt.b%rampt1))then
        ramptw = 0.5d0*(1d0 - dcos(pi*(t - b%rampt0)/b%rampdt) )
      endif

      call wv%getEta(t,0d0,0d0,eta)
      call wv%getPQ(t,0d0,0d0,eta,p,q)

      b%data(i+1,1)=t
      b%data(i+1,2)=eta * ramptw
      b%data(i+1,3)=p * ramptw
      b%data(i+1,4)=q * ramptw
      !write(201,'(4F15.6)')b%data(i+1,1:4)
    enddo    

    b%datadtt(:,1) = b%data(:,1)
    call setCubicSpline(b%numP, 4, b%data, b%datadtt)
    b%posI=1  

  end subroutine initAiryFile
!!----------------------End initAiryFile-----------------------!!



!!-----------------------initStokes2File-----------------------!!
  subroutine initStokes2File(b,dt,inTotT,inT,inD,inH,inAng, &
    rampt0, rampt1)
  use stokes2WaveModule
  implicit none

    class(wvFileType),intent(inout)::b
    integer(kind=C_K1)::i
    real(kind=C_K2),intent(in)::inTotT,inT,inD,inH,dt,inAng
    real(kind=C_K2),intent(in)::rampt0, rampt1
    real(kind=C_K2)::t, eta, p, q, totT, ramptw
    type(stokes2Type)::wv    

    ! stokes2Type(T,d,H,X0,Y0,thDeg)
    wv=stokes2Type(inT,inD,inH,0d0,0d0,inAng)

    write(9,'(" [INF] ",3A15)')'T','L','d'
    write(9,'(" [---] ",3F15.6)')wv%T,wv%L,wv%d
    write(9,'(" [INF] ",A15)')'kh'
    write(9,'(" [---] ",F15.6)')wv%k*wv%d
    write(9,'(" [INF] At 0.00 WavePeriods")')
    write(9,'(" [---] ",3A15)')'Eta','P','Q'
    call wv%getEta(0d0,wv%x0,wv%y0,eta)
    call wv%getPQ(0d0,wv%x0,wv%y0,eta,p,q)
    write(9,'(" [---] ",3F15.6)')eta,p,q
    write(9,'(" [INF] At 2.25 WavePeriods")')
    write(9,'(" [---] ",3A15)')'Eta','P','Q'
    call wv%getEta(2.25d0*wv%T,wv%x0,wv%y0,eta)
    call wv%getPQ(2.25d0*wv%T,wv%x0,wv%y0,eta,p,q)
    write(9,'(" [---] ",3F15.6)')eta,p,q

    totT=1.2d0*inTotT
    b%numP=floor(totT/dt)+2
    b%rampt0 = rampt0
    b%rampt1 = rampt1
    b%rampdt = (rampt1 - rampt0)

    allocate(b%data(b%numP,4), b%datadtt(b%numP,4))
    do i=0,b%numP-1
      t=i*dt

      !Time Ramp
      ramptw = 1d0
      if(t.le.b%rampt0)then
        ramptw = 0d0
      else if(t.ge.b%rampt1)then
        ramptw = 1d0
      else if((t.gt.b%rampt0) .and. (t.lt.b%rampt1))then
        ramptw = 0.5d0*(1d0 - dcos(pi*(t - b%rampt0)/b%rampdt) )
      endif

      call wv%getEta(t,0d0,0d0,eta)
      call wv%getPQ(t,0d0,0d0,eta,p,q)

      b%data(i+1,1)=t
      b%data(i+1,2)=eta * ramptw
      b%data(i+1,3)=p * ramptw
      b%data(i+1,4)=q * ramptw
      !write(201,'(4F15.6)')b%data(i+1,1:4)
    enddo    

    b%datadtt(:,1) = b%data(:,1)
    call setCubicSpline(b%numP, 4, b%data, b%datadtt)
    b%posI=1  

  end subroutine initStokes2File
!!---------------------End initStokes2File---------------------!!




!!---------------------------getEta----------------------------!!
  ! subroutine getEta(b,rTime,eta)
  ! implicit none

  !   class(wvFileType),intent(in)::b

  !   integer(kind=C_K1)::k,posI

  !   real(kind=C_K2),intent(in)::rTime
  !   real(kind=C_K2),intent(out)::eta

  !   posI=2
    
  !   if((rTime.gt.b%data(b%numP,1)).or.&
  !     (rTime.lt.b%data(1,1)))then
  !     write(9,*)"[ERR] Wave query exceeds supplied time range"
  !     write(9,'( "[---] ",F15.6)')rTime
  !     write(9,'( "[---] ",2F15.6)')b%data(1,1),b%data(b%numP,1)
  !     stop
  !   endif

  !   do while(b%data(posI+1,1).le.rTime)
  !     posI=posI+1
  !   enddo

  !   ! do while(b%data(posI-1,1).gt.rTime)
  !   !   if(posI.eq.2)exit
  !   !   posI=posI-1
  !   ! enddo

  !   k=2
  !   eta=b%data(posI,k) &
  !     + (b%data(posI+1,k)-b%data(posI-1,k)) &
  !     / (b%data(posI+1,1)-b%data(posI-1,1)) &
  !     * (rTime-b%data(posI,1))

  !   !write(222,'(3F20.6)')rTime,b%data(posI,1),rTime-b%data(posI,1)

  ! end subroutine getEta



  subroutine getEta(b,rTime,eta)
  implicit none

    class(wvFileType),intent(in)::b

    integer(kind=C_K1)::k, posI, i

    real(kind=C_K2),intent(in)::rTime
    real(kind=C_K2),intent(out)::eta
    real(kind=C_K2)::dx, cA, cB, cC, cD

    ! if((rTime.ge.b%data(b%posI,1)) .and. &
    !   (rTime.lt.b%data(b%posI+1,1)))then
      
    !   posI=b%posI      

    ! elseif((rTime.ge.b%data(b%posI+1,1)) .and. &
    !   (rTime.lt.b%data(b%posI+2,1)))then
      
    !   b%posI = b%posI + 1
    !   posI=b%posI

    ! elseif(b%posI.ge.2)then
    !   if((rTime.ge.b%data(b%posI-1,1)) .and. &
    !     (rTime.lt.b%data(b%posI,1)))then

    !     b%posI = b%posI - 1
    !     posI=b%posI
    !   endif
    
    ! else
      
    !   do i = 1, b%numP
    !     if(b%data(i,1).gt.rTime) exit      
    !   enddo
    !   posI=i-1      
    !   if((posI.eq.0) .or. (rTime.gt.b%data(b%numP,1)) )then        
    !     write(9,*)"[ERR] Wave query exceeds supplied time range"
    !     write(9,'( "[---] ",F15.6)')rTime
    !     write(9,'( "[---] ",2F15.6)')b%data(1,1),b%data(b%numP,1)
    !     stop
    !   endif
    !   b%posI=posI
    !   write(9,*)"[INF] Looped for waveInput"
    ! endif  

    do i = 1, b%numP
      if(b%data(i,1).gt.rTime) exit      
    enddo
    posI=i-1      
    if((posI.eq.0) .or. (rTime.gt.b%data(b%numP,1)) )then        
      write(9,*)"[ERR] Wave query exceeds supplied time range"
      write(9,'( "[---] ",F15.6)')rTime
      write(9,'( "[---] ",2F15.6)')b%data(1,1),b%data(b%numP,1)
      stop
    endif  

    dx = (b%data(posI+1,1) - b%data(posI,1))
    cA = (b%data(posI+1,1) - rTime) / dx
    cB = 1d0 - cA
    cC = (cA**3 - cA)*dx*dx/6d0
    cD = (cB**3 - cB)*dx*dx/6d0          

    k=2
    eta = (cA * b%data(posI,k)) + (cB * b%data(posI+1,k)) &
      + (cC * b%datadtt(posI,k)) + (cD * b%datadtt(posI+1,k))    

    !write(222,'(3F20.6)')rTime,b%data(posI,1),rTime-b%data(posI,1)

  end subroutine getEta
!!-------------------------End getEta--------------------------!!



!!----------------------------getPQ----------------------------!!
  subroutine getPQ(b,rTime,p,q)
  implicit none

    class(wvFileType),intent(in)::b

    integer(kind=C_K1)::k, posI, i

    real(kind=C_K2),intent(in)::rTime
    real(kind=C_K2),intent(out)::p,q
    real(kind=C_K2)::dx, cA, cB, cC, cD

    do i = 1, b%numP
      if(b%data(i,1).gt.rTime) exit      
    enddo
    posI=i-1      
    if((posI.eq.0) .or. (rTime.gt.b%data(b%numP,1)) )then        
      write(9,*)"[ERR] Wave query exceeds supplied time range"
      write(9,'( "[---] ",F15.6)')rTime
      write(9,'( "[---] ",2F15.6)')b%data(1,1),b%data(b%numP,1)
      stop
    endif 

    dx = (b%data(posI+1,1) - b%data(posI,1))
    cA = (b%data(posI+1,1) - rTime) / dx
    cB = 1d0 - cA
    cC = (cA**3 - cA)*dx*dx/6d0
    cD = (cB**3 - cB)*dx*dx/6d0          

    k=3
    p = (cA * b%data(posI,k)) + (cB * b%data(posI+1,k)) &
      + (cC * b%datadtt(posI,k)) + (cD * b%datadtt(posI+1,k)) 

    k=4
    q = (cA * b%data(posI,k)) + (cB * b%data(posI+1,k)) &
      + (cC * b%datadtt(posI,k)) + (cD * b%datadtt(posI+1,k)) 
    
  end subroutine getPQ
!!--------------------------End getPQ--------------------------!!



!!----------------------chkWaveFileInput-----------------------!!
  subroutine chkWaveFileInput(b,rT0,rT1,dt)
  implicit none

    class(wvFileType),intent(in)::b

    integer(kind=C_K1)::i,j,mf,nT

    real(kind=C_K2),intent(in)::rT0,rT1,dt
    real(kind=C_K2)::eta,p,q,t

    open(newunit=mf,file='chkWaveFileInput.dat')

    nT=floor((rT1-rT0)/dt)+1
    do i=0,nT
      t=rT0+i*dt
      call b%getEta(t,eta)
      call b%getPQ(t,p,q)
      write(mf,'(4F20.8)')t,eta,p,q
    enddo
    close(mf)


  end subroutine chkWaveFileInput
!!--------------------End chkWaveFileInput---------------------!!



!!-----------------------setCubicSpline------------------------!!
  subroutine setCubicSpline(nn, nc, da, dadtt)
  implicit none
    
    integer(kind=C_K1),intent(in)::nn, nc
    integer(kind=C_K1)::i,k    
    real(kind=C_K2),intent(in)::da(nn, nc)
    real(kind=C_K2),intent(out)::dadtt(nn, nc)
    real(kind=C_K2),allocatable::mA(:,:),rhs(:,:)
    

    ! mA (n x n)
    ! [b1 c1 ----------------------]
    ! [a2 b2 c2 -------------------]
    ! [-- a3 b3 c3 ----------------]
    ! [----- a4 b4 c4 -------------]
    ! [----------------------------]
    ! [---------------------- an bn]
    !
    ! Convert to
    ! 
    ! mA (n x 3)
    ! [ 0 b1 c1]
    ! [a2 b2 c2]
    ! [a3 b3 c3]
    ! [a4 b4 c4]
    ! [--------]
    ! [an bn 0 ]

    ! rhs
    ! corresponding to 
    ! [x, y, th]

    write(9,*)
    write(9,'(" [INF] waveInput Cublic Spine Coeffs")')    
    
    allocate(mA(nn,3), rhs(nn,nc-1))
    mA=0d0
    rhs=0d0

    ! natural spline setup
    mA(1,2)=1d0 !b1
    mA(nn,2)=1d0 !bn
    do i=2,nn-1
      mA(i,1) = (da(i,1)-da(i-1,1))/6d0
      mA(i,2) = (da(i+1,1)-da(i-1,1))/3d0
      mA(i,3) = (da(i+1,1)-da(i,1))/6d0

      do k = 2,nc
        rhs(i,k-1) = (da(i+1,k) - da(i,k)) &
          / (da(i+1,1) - da(i,1)) - &
          (da(i,k) - da(i-1,k)) &
          / (da(i,1) - da(i-1,1)) 
      enddo
      
    enddo    

    do k = 2,nc
      call triDiagSolver(nn, mA, rhs(:,k-1), dadtt(:,k))
    enddo

    deallocate(mA, rhs)
    write(9,'(" [---] Done")')    
    write(9,*)
    
  end subroutine setCubicSpline
!!---------------------End setCubicSpline----------------------!!



!!------------------------triDiagSolver------------------------!!
  subroutine triDiagSolver(n,A,B,x)
  implicit none

    integer(kind=C_K1),intent(in)::n
    real(kind=C_K2),intent(in)::A(n,3),B(n)
    real(kind=C_K2),intent(out)::x(n)

    integer(kind=C_K1)::i    
    real(kind=C_K2),allocatable::prC(:),prD(:)    

    ! A (n x n)
    ! [b1 c1 ----------------------]
    ! [a2 b2 c2 -------------------]
    ! [-- a3 b3 c3 ----------------]
    ! [----- a4 b4 c4 -------------]
    ! [----------------------------]
    ! [---------------------- an bn]
    !
    ! Converted to following in input itself
    ! 
    ! A (n x 3)
    ! [ 0 b1 c1]
    ! [a2 b2 c2]
    ! [a3 b3 c3]
    ! [a4 b4 c4]
    ! [--------]
    ! [an bn 0 ]

    ! [ai bi ci][xi] = [di]      

    allocate(prC(n), prD(n))

    prC(1) = A(1,3)/A(1,2)
    prD(1) = B(1)/A(1,2)

    do i = 2,n-1
      prC(i) = A(i,3) / ( A(i,2) - A(i,1)*prC(i-1) )      
    enddo

    do i = 2,n      
      prD(i) = ( B(i) - A(i,1)*prD(i-1) ) &
        / ( A(i,2) - A(i,1)*prC(i-1) )
    enddo

    x(n) = prD(n)

    do i = n-1,1,-1
      x(i) = prD(i) - prC(i)*x(i+1)
    enddo

    deallocate(prC, prD)

  end subroutine triDiagSolver
!!----------------------End triDiagSolver----------------------!!


end module waveFileModule
!!-----------------------End waveFileModule------------------------!!