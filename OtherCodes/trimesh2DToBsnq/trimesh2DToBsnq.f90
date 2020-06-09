!! Edit 
!! 2018-12-31
!! 1430
program trimesh2DToBsnq
implicit none

  integer(kind=4)::i,j,k,l,i2,j2,k2,l2,fl1,fl2,fl3
  integer(kind=4)::np,nele,nbndty,nedg,nedgMax,nbnd
  integer(kind=8)::timeVal(8)  
  logical::ex
  integer(kind=4),allocatable::conn(:,:),bndp(:,:),edg(:,:)
  integer(kind=4),allocatable::bnd(:,:)
  real(kind=8),allocatable::coor(:,:)
  character(len=200)::probname,text,outname

  call getarg(1,probname)  
  if(len_trim(probname).lt.1) then
    write(*,'(A)',advance='no')"Enter Problem Name:"
    read(*,*)probname
  endif
  write(*,201)"[INF] ","Problem Name: ",trim(probname)  

  fl1=10
  inquire(file=trim(probname),exist=ex)
  if(ex)then
    open(fl1,file=trim(probname))
  else
    write(*,'(a7,a)')'[ERR] ','Missing mesh file'
  endif

!!--------------------------Reading---------------------------!!
  read(fl1,*)np
  allocate(coor(np,2))

  do i=1,np
    read(fl1,*)l2,coor(i,1:2)
  enddo

  read(fl1,*)nele
  allocate(conn(nele,3))
  do i=1,nele
    read(fl1,*)l2,conn(i,1:3)
  enddo

  read(fl1,*)
  read(fl1,*)nbndty
  nbndty=nbndty+1
  allocate(bndp(nbndty,-1:np))
  bndp=0

  !! Named Boundaries
  do k=1,nbndty-1
    read(fl1,*)
    read(fl1,*)bndp(k,-1)
    read(fl1,*)
    read(fl1,*)bndp(k,0)
    read(fl1,*)
    l2=floor(1d0*bndp(k,0)/10d0)
    j=0
    do i=1,l2
      read(fl1,*)bndp(k,j+1:j+10)
      j=j+10
    enddo
    read(fl1,*)bndp(k,j+1:bndp(k,0))
  enddo

  !! Remaining boundaries
  read(fl1,*)nbnd
  k=nbndty
  bndp(k,-1)=13
  i2=0
  do i=1,nbnd
    read(fl1,*)l2,k2
    j2=0
    do j=1,nbndty-1
      j2=j2+count(bndp(j,1:bndp(j,0)).eq.k2)
    enddo
    if(j2.eq.0)then
      i2=i2+1
      bndp(k,i2)=k2
    endif
  enddo
  bndp(k,0)=i2

  write(*,'(a7,2a15)')"[INF] ","BndTyp","NumBndTyp"
  do i=1,nbndty
    write(*,'(a7,2i15)')'[---] ',bndp(i,-1:0)
  enddo
  write(*,'(a7,a15,i15)')'[---] ','Sum Total',sum(bndp(:,0))
  write(*,'(a7,a15,i15)')'[---] ','NumBnd Read',nbnd
  if(nbnd.ne.sum(bndp(:,0)))then
    write(*,'(a7,a)')'[ERR] ','Check numBnd and bndTyp Loc1'
    write(*,'(a7,a,i15)')'[---] ','Num Bnd',nbnd
    stop
  endif

  nedgMax=3*nele
  allocate(edg(nedgMax,5))
  !! edg(i,1:5) = (/ lowerNode, higherNode, mul, ele1, ele2 /)
  !! lower node   : the node with lesser node number
  !! higher node  : the node with higher node number
  !! mul          : anticlk = 1, clk = 1 only imp for bndEdges    
  call findEdg(nele,nedg,nedgMax,conn,edg)
  nbnd=count(edg(:,3).ne.0)
  write(*,'(a7,4a15)')'[INF] ','Elements','Nodes','Edges','Bnd'
  write(*,'(a7,4i15)')'[---] ',nele,np,nedg,nbnd
  !! NumBndEdges = NumBndPoi, because we have closed boundaries
  if(nbnd.ne.sum(bndp(:,0)))then
    write(*,'(a7,a)')'[ERR] ','Check numBnd and bndTyp Loc2'
    write(*,'(a7,a,i15)')'[---] ','Num Bnd',nbnd
    stop
  endif

  allocate(bnd(nbnd,4))
  k=0
  bnd=0
  do i=1,nedg
    if(edg(i,3).eq.0)cycle
    k=k+1
    if(edg(i,3).eq.1)then
      bnd(k,:)=(/ edg(i,1),edg(i,2),edg(i,4),0 /)
    else
      bnd(k,:)=(/ edg(i,2),edg(i,1),edg(i,4),0 /)
    endif
  enddo
  if(k.ne.nbnd)then
    write(*,'(a7,a,2i15)')'[ERR] ','Check nbnd.ne.k',nbnd,k
    stop
  endif
  

  do i=1,nbnd
    k=bnd(i,1)
    k2=bnd(i,2)

    l=0
    l2=0
    do j=1,nbndty
      if(count(bndp(j,1:bndp(j,0)).eq.k).ne.0)then
        if(l.ne.0)then
          write(*,'(a7,a)')'[ERR] ','Bnd node in multiple types'
          write(*,'(a7,i15)')'[---] ',k
          stop
        endif
        l=bndp(j,-1)
      endif
      if(count(bndp(j,1:bndp(j,0)).eq.k2).ne.0)then
        if(l2.ne.0)then
          write(*,'(a7,a)')'[ERR] ','Bnd node in multiple types'
          write(*,'(a7,i15)')'[---] ',k2
          stop
        endif
        l2=bndp(j,-1)
      endif
    enddo
    if((l*l2).eq.0)then
      write(*,'(a7,a)')'[ERR] ','BndEdge nodes not found in bndp'
      write(*,'(a7,3a15)')'[---] ','BndNo','Node1','Node2'
      write(*,'(a7,3i15)')'[---] ',i,k,k2
      stop
    endif
    bnd(i,4)=max(l,l2)
    !! Done assuming that + are inlets(11) and - are walls(13)
    !!    +------------- 
    !!    +
    !!    +
    !!    +
    !!    +-------------
  enddo
  i=count(bnd(:,4).eq.0)
  if(i.ne.0)then
    write(*,'(a7,a,i15,a)')'[ERR] ','BndTyp not defined for ',&
      i,'edges'
  endif
!!------------------------End Reading-------------------------!!

!!--------------------------Writing---------------------------!!
  call date_and_time(VALUES=timeVal)
  write(outname,'("Mesh_",I4.4,4I2.2,".plt")')timeVal(1:3),timeVal(5:6)
  !write(outname,'("Mesh_",I4.4,4I2.2,".plt")')timeVal(1:3),timeVal(5),2
  write(*,201)"[INF] ",'Output Name: ',trim(outname)
  fl2=11
  open(fl2,file=trim(outname))

  !! Nodes only
  write(outname,'("Mesh_",I4.4,4I2.2,"_nodes.dat")')&
    timeVal(1:3),timeVal(5:6)
  fl3=12
  open(fl3,file=trim(outname))
  write(fl3,'("Number of Nodes")')
  write(fl3,'(i15)')np
  write(fl3,'("Nodes")')
  do i=1,np
    write(fl3,'(2f20.6)')coor(i,:)
  enddo

  !! Mesh file
  write(fl2,'(3a15)')'Elements','Nodes','Edges'
  write(fl2,'(3i15)')nele,np,nedg
  write(fl2,'(2a15)')'BndSides','BndTypes'
  write(fl2,'(2i15)')nbnd,nbndty
  write(fl2,'("Nodes")')
  do i=1,np
    write(fl2,'(2f20.6)')coor(i,:)
  enddo  
  write(fl2,'("Elements")')
  do i=1,nele
    write(fl2,'(3i15)')conn(i,:)
  enddo
  do k=1,nbndty
    k2=bndp(k,-1)
    write(fl2,'("Boundary")')
    l2=count(bnd(:,4).eq.k2)
    write(fl2,'(2i15)')k2,l2
    do j=1,nbnd
      if(bnd(j,4).eq.k2)then 
        write(fl2,'(3i15)')bnd(j,1:3)
      endif
    enddo
  enddo
  write(fl2,'("Depth")')

  
!!------------------------End Writing-------------------------!!


  close(fl1)
  close(fl2)
  close(fl3)
  201 format(a7,a20,a)

end program trimesh2DToBsnq

subroutine findEdg(nele,nedg,nedgMax,conn,edg)
implicit none

  integer(kind=4)::i,j,k,l,i2,j2,k2,iel,p1,p2,lmul  
  integer(kind=4)::n(3),perc1,perc2
  integer(kind=4),intent(in)::nele,nedgMax,conn(nele,3)
  integer(kind=4),intent(out)::nedg,edg(nedgMax,5)  

  nedg=0
  edg=0

  write(*,'(a7,a)')'[MSG] ','Entering findEdg'

  perc1=0
  do iel=1,nele
    perc2=floor(10d0*iel/nele)
    if(perc1.ne.perc2)then
      write(*,'(a7,i15,"%")')'[---] ',perc2*10
      perc1=perc2
    endif

    n=conn(iel,:)

    do i=1,3
      i2=mod(i,3)+1
      if(n(i).lt.n(i2))then
        lmul=1
        p1=n(i)
        p2=n(i2)
      else
        lmul=-1
        p1=n(i2)
        p2=n(i)
      endif

      k=count(edg(1:nedg,1).eq.p1)
      if(k.eq.0)then

        nedg=nedg+1
        edg(nedg,:)=(/ p1,p2,lmul,iel,0 /)

      else

        k2=0
        do j=1,nedg
          if(edg(j,1).eq.p1)then
            if(edg(j,2).eq.p2)then
              k2=1
              edg(j,3)=0
              edg(j,5)=iel
              exit
            endif
          endif
        enddo
        if(k2.eq.0)then
          nedg=nedg+1
          edg(nedg,:)=(/ p1,p2,lmul,iel,0 /)
        endif

      endif
    enddo

  enddo

  write(*,'(a7,a,i15)')'[INF] ','Num Edge: ',nedg


end subroutine