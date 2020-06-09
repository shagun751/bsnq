include 'output.f90'

program genQuadMesh
implicit none

  integer(kind=4)::i,j,k,l,i2,j2,k2,nx,ny,en(4)
  integer(kind=4)::npoi,nele,nbnd,nedg,bndt1
  integer(kind=4),allocatable::conn(:,:),mabnd(:,:)
  real(kind=8)::domx(2),domy(2)
  real(kind=8)::tmpr1,tmpr2,tmpr3,tmpr4,tmpr5,dx,dy
  real(kind=8),allocatable::corx(:),cory(:)
  real(kind=8),allocatable::dep(:)
  character(len=100)::probname

  probname='test'

  domx(1)=79d0
  domy(1)=10d0
  domx(2)=83d0
  domy(2)=14d0
  dx=0.05d0
  dy=0.05d0
  nx=int((domx(2)-domx(1))/dx,4)
  ny=int((domy(2)-domy(1))/dy,4)
  write(*,*)nx,ny

  npoi=(nx+1)*(ny+1)
  nele=nx*ny
  nedg=2*nx*ny+nx+ny
  nbnd=2*(nx+ny)

  allocate(corx(npoi),cory(npoi),conn(nele,4))
  allocate(mabnd(nbnd,3),dep(npoi))

  !! Nodes
  k=0
  do i=1,nx+1
    do j=1,ny+1
      k=k+1
      corx(k)=(i-1)*dx+domx(1)
      cory(k)=(j-1)*dy+domy(1)
    enddo
  enddo
  write(*,*)k,npoi
  ! do i=1,npoi
  !   write(*,'(2F15.6)')corx(i),cory(i)
  ! enddo

  !! Elements
  k=0 
  do i=1,nx
    do j=1,ny
      k=k+1
      k2=(i-1)*(ny+1)+j
      en(1)=k2
      en(2)=k2+ny+1
      en(3)=k2+ny+2
      en(4)=k2+1
      conn(k,:)=en
    enddo
  enddo
  write(*,*)k,nele
  ! do i=1,nele
  !   write(*,'(4I15)')conn(i,:)
  ! enddo

  !! Boundaries
  k=0
  do i=1,ny
    k=k+1
    mabnd(k,1)=i+1
    mabnd(k,2)=i
    mabnd(k,3)=i
  enddo
  bndt1=k
  do i=1,nx
    k=k+1
    j=(i-1)*ny+1
    mabnd(k,1)=conn(j,1)
    mabnd(k,2)=conn(j,2)
    mabnd(k,3)=j
  enddo
  do i=1,ny
    k=k+1
    j=nx*ny-ny+i
    mabnd(k,1)=conn(j,2)
    mabnd(k,2)=conn(j,3)
    mabnd(k,3)=j
  enddo
  do i=nx,1,-1
    k=k+1
    j=i*ny
    mabnd(k,1)=conn(j,3)
    mabnd(k,2)=conn(j,4)
    mabnd(k,3)=j
  enddo
  write(*,*)k,nbnd

  
  ! dep=100d0
  ! do i=1,nbnd
  !   dep(mabnd(i,1))=12
  !   dep(mabnd(i,2))=12
  ! enddo

  ! tmpr5=0.12d0  !! Continental slop generally
  ! tmpr4=220850
  ! dep=3500d0
  ! do i=1,npoi
  !   tmpr1=abs(corx(i)-250000)
  !   if(tmpr1.gt.tmpr4)then
  !     tmpr1=tmpr1-tmpr4
  !     dep(i)=3500d0-0.12d0*tmpr1
  !   endif
  ! enddo

  dep=500
  
  call out4NXML(probname,npoi,npoi,nele,21,0,&
  conn,corx,cory,dep,dep,dep,dep)

  !! Mesh File
  write(11,'(2A15)')"Elements","Nodes"
  write(11,'(2I15)')nele,npoi
  write(11,'(2A15)')"BndSides","BndTypes"
  write(11,'(2I15)')nbnd,2
  write(11,'("Nodes")')
  do i=1,npoi
    write(11,'(2F20.8)')corx(i),cory(i)
  enddo
  write(11,'("Elements")')
  do i=1,nele
    write(11,'(4I15)')conn(i,:)
  enddo
  write(11,'("Boundary")')
  write(11,'(2I15)')12,bndt1
  do i=1,bndt1
    write(11,'(3I15)')mabnd(i,:)
  enddo
  write(11,'("Boundary")')
  write(11,'(2I15)')12,nbnd-bndt1
  do i=bndt1+1,nbnd
    write(11,'(3I15)')mabnd(i,:)
  enddo
  write(11,'("Depth")')
  do i=1,npoi
    write(11,'(F15.6)')dep(i)
  enddo
  write(11,'("####")')


end program genQuadMesh