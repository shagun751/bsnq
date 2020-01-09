module shipInfo
use basicVars
implicit none

  real(kind=C_K2)::cb,cl,al
  real(kind=C_K2)::L,B,P0
  real(kind=C_K2)::cx,cy,v,th
  real(kind=C_K2),allocatable::posData(:,:)

  integer(kind=C_K1)::posN,posI

end module shipInfo


subroutine shipInput(probname,rtime)
use shipInfo
implicit none
    
  real(kind=C_K2),intent(in)::rtime
  character(len=100),intent(in)::probname
  character(len=100)::text
  logical::ex

  text=trim(probname)//".pos"
  inquire(file=trim(text),exist=ex)
  if(.not.ex)then
    write(9,'(a7,a)')'[ERR] ','No ship position file of type .pos'
    write(9,'(a7,a,a)')'[---] ','Missing ',trim(text)
    stop
  endif
  open(5,file=trim(text))

  read(5,*,end=11,err=11)text
  read(5,*,end=11,err=11)text
  read(5,*,end=11,err=11)cl,cb,al
  read(5,*,end=11,err=11)text
  read(5,*,end=11,err=11)L,B,P0
  read(5,*,end=11,err=11)text
  read(5,*,end=11,err=11)posN
  read(5,*,end=11,err=11)text
  read(5,*,end=11,err=11)cx,cy
  read(5,*,end=11,err=11)text
  allocate(posData(posN,3))
  posI=0
  do i=1,posN
    read(5,*,end=11,err=11)posData(i,:)
    if(rtime.ge.posData(i,1))posI=i
  enddo
  if(posI.eq.0)then
    write(9,'(a7,a,f15.6,a)')'[ERR] ',&
      'Ship position data start at',&
      posData(1,1),' s'
    write(9,'(a7,a,f15.6,a)')'[ERR] ',&
      'Simulation starts at',rtime,' s'
    stop
  endif

  posData(:,3)=posData(:,3)*pi/180
  P0=P0*grav*rhoW

  write(9,*)
  write(9,'(a7,a)')'[INF] ','Ship shape info'
  write(9,'(a7,3a15)')'[---] ','cl','cb','al'
  write(9,'(a7,3f15.6)')'[---] ',cl,cb,al
  write(9,'(a7,3a15)')'[---] ','L','B','P0'
  write(9,'(a7,3f15.6)')'[---] ',L,B,P0 
  write(9,'(a7,a)')'[INF] ','Ship position'
  write(9,'(a7,3a15)')'[---] ','T(s)','V(m/s)','th(rad)'
  do i=1,posN
    write(9,'(a7,3f15.6)')'[---] ',posData(i,:)
  enddo

  goto 15
  11 write(9,'(a7,a)')'[ERR] ','Check ship position file format'
  stop
  15 close(5)

  write(9,*)
end subroutine shipInput


subroutine shipPosition(rtime,dt)
use shipInfo
implicit none
  
  real(kind=C_K2),intent(in)::rtime,dt
  real(kind=C_K2)::val0,val1,t0,t1,tVal
  real(kind=C_K2)::vx,vy
  
  if(rtime.gt.posData(posN,1))then
    posI=posN
    v=posData(posN,2)
    th=posData(posN,3)
    return
  endif
  if(rtime.gt.posData(posI+1,1)) posI=posI+1

  t0=posData(posI,1)
  t1=posData(posI+1,1)
  tVal=rtime-dt/2d0
  
  val0=posData(posI,2)
  val1=posData(posI+1,2)
  v=(val1-val0)/(t1-t0)*(tVal-t0)+val0

  val0=posData(posI,3)
  val1=posData(posI+1,3)
  th=(val1-val0)/(t1-t0)*(tVal-t0)+val0

  vx=v*dcos(th)
  vy=v*dsin(th)
  cx=cx+vx*dt
  cy=cy+vy*dt

  write(9,201)'SHPP',cx,cy
  write(9,201)'SHPV',v,(th*180d0/pi)

  201 format('      |',a5,2f15.6)
end subroutine shipPosition


subroutine shipPress(npt,coord,presDiv)
use shipInfo
implicit none

  integer(kind=C_K1),intent(in)::npt
  real(kind=C_K2),intent(in) ::coord(npt,2)
  real(kind=C_K2),intent(out)::presDiv(npt,2)    
  real(kind=C_K2)::lx,ly,ly2,rRef
  real(kind=C_K2)::cost,sint,x,y  

  presDiv=0d0

  rRef=(L/2d0*1.2d0)**2
  cost=dcos(th)
  sint=dsin(th)
  
  do i=1,npt
    x=coord(i,1)-cx
    y=coord(i,2)-cy
    tmpr3=(x**2 + y**2)
    if(tmpr3.le.rRef)then
      lx=(+x*cost + y*sint)/L
      ly=(-x*sint + y*cost)/B
      ! lx=(-x*cost - y*sint)/L
      ! ly=(+x*sint - y*cost)/B
      if(abs(lx).gt.0.5d0)cycle
      if(abs(ly).gt.0.5d0)cycle
      ly2=ly**2
      tmpr4=1d0/dexp(al*(ly2))
      tmpr1=(-4d0*cl*P0*(lx**3)&
        *(1d0 - cb*(ly2)))*tmpr4
      tmpr2=(-2d0*P0)*tmpr4 &
        *(-1d0 + cl*(lx**4))*ly &
        *(-cb + al*(-1d0 + cb*ly2))
      tmpr1=tmpr1/L
      tmpr2=tmpr2/B
      presDiv(i,1)=tmpr1*cost - tmpr2*sint
      presDiv(i,2)=tmpr1*sint + tmpr2*cost
      ! presDiv(i,1)=-tmpr1*cost + tmpr2*sint
      ! presDiv(i,2)=-tmpr1*sint - tmpr2*cost
    endif
  enddo

end subroutine shipPress

subroutine shipInitEta(npl,npt,coord,eta,p,q)
use shipInfo
implicit none

  integer(kind=C_K1),intent(in)::npl,npt
  real(kind=C_K2),intent(in) ::coord(npt,2)
  real(kind=C_K2),intent(out)::eta(npl),p(npt),q(npt)  
  real(kind=C_K2)::lx,ly,ly2,rRef
  real(kind=C_K2)::cost,sint,x,y  
  

  p=0d0
  q=0d0
  eta=0d0

  rRef=(L/2d0*1.2d0)**2
  cost=dcos(th)
  sint=dsin(th)
  
  do i=1,npt
    x=coord(i,1)-cx
    y=coord(i,2)-cy
    tmpr3=(x**2 + y**2)
    if(tmpr3.le.rRef)then
      lx=(x*cost + y*sint)/L
      ly=(-x*sint + y*cost)/B
      if(abs(lx).gt.0.5d0)cycle
      if(abs(ly).gt.0.5d0)cycle
      ly2=ly**2
      eta(i)=-0.5d0*(1d0-cl*(lx**4)) &
        *(1d0-cb*ly2)*dexp(-al*ly2)
    endif
  enddo


end subroutine shipInitEta