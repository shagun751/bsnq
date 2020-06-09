program quadmesh2DToBsnq
implicit none

  integer(kind=4)::i,j,k,l,i2,j2,k2,l2,fl1,fl2,fl3,fli
  integer(kind=4)::np,nele,nbndty,nbnd,nbnd2
  integer(kind=4)::p(4),pSt(4),eleBnd(4,4),nEleBnd
  integer(kind=8)::timeVal(8)  
  logical::ex
  integer(kind=4),allocatable::conn(:,:),bndp(:,:)
  integer(kind=4),allocatable::bnd(:,:)
  real(kind=8),allocatable::coor(:,:)
  real(kind=8)::sc1cx,sc1cy,sc1lon,sc1lat,sc2cx,sc2cy,sc2lon,sc2lat
  real(kind=8)::scdx,scdy
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

!!---------------------------Input----------------------------!!
  fli=13
  open(fli,file='Scale.txt')
  read(fli,*)text
  read(fli,*)text
  read(fli,*)sc1cx,sc1cy
  read(fli,*)text
  read(fli,*)sc1lon,sc1lat
  read(fli,*)text
  read(fli,*)text
  read(fli,*)sc2cx,sc2cy
  read(fli,*)text
  read(fli,*)sc2lon,sc2lat  
  scdx=(sc2lon-sc1lon)/(sc2cx-sc1cx)
  scdy=(sc2lat-sc1lat)/(sc2cy-sc1cy)
  write(*,*)
  write(*,'(a7,a)')'[INF] ','Scaling'
  write(*,'(a7,a)')'[---] ','Bottom Left'
  write(*,'(a7,a)')'[---] ','Mesh Cx Cy'
  write(*,'(a7,2F15.6)')'[---] ',sc1cx,sc1cy
  write(*,'(a7,a)')'[---] ','Map Lon Lat'
  write(*,'(a7,2F15.6)')'[---] ',sc1lon,sc1lat
  write(*,'(a7,a)')'[---] ','Top Right'
  write(*,'(a7,a)')'[---] ','Mesh Cx Cy'
  write(*,'(a7,2F15.6)')'[---] ',sc2cx,sc2cy
  write(*,'(a7,a)')'[---] ','Map Lon Lat'
  write(*,'(a7,2F15.6)')'[---] ',sc2lon,sc2lat
  write(*,'(a7,a)')'[---] ','scDx, scDy'
  write(*,'(a7,2F15.6)')'[---] ',scdx,scdy
  write(*,*)
!!-------------------------End Input--------------------------!!

!!--------------------------Reading---------------------------!!
  read(fl1,*)np
  allocate(coor(np,2))

  do i=1,np
    read(fl1,*)l2,coor(i,1:2)
  enddo

  read(fl1,*)nele
  allocate(conn(nele,4))
  do i=1,nele
    read(fl1,*)l2,conn(i,1:4)
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

  write(outname,'("temp")')
  call out4NXML(outname,np,np,nele,5,0,conn,&
    coor(:,1),coor(:,2),coor(:,1),coor(:,2),coor(:,1),coor(:,2))

  !! Bnd detect
  nbnd2=0
  allocate(bnd(nbnd,4))
  do i=1,nele
    p=conn(i,:)
    nEleBnd=0
    eleBnd=0
    
    do i2=0,3
      k=p(i2+1)
      k2=p(mod(i2+1,4)+1)

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
      ! if((k.eq.918).and.(k2.eq.915)&
      !   .or.(k.eq.915).and.(k2.eq.918)) cycle
      !if((k.eq.5300).or.(k2.eq.5300))cycle

      if((l*l2).ne.0)then
        nEleBnd=nEleBnd+1
        eleBnd(nEleBnd,:)=(/ k,k2,i,max(l,l2) /)
      endif
      !! Done assuming that + are inlets(11) and - are walls(13)
      !!    +------------- 
      !!    +
      !!    +
      !!    +
      !!    +-------------      
    enddo

    if(nEleBnd.eq.4)then
      write(*,'(a7,a,i10)')'[ERR] ','Note 4 bndNodes in ele ',i
      write(*,'(a7,a)')'[---] ','The algortihm will fail'
      write(*,'(a7,a)')'[---] ','Location'
      do j=1,4
        write(*,'(a7,I10,2F15.6)')'[---] ',p(j),coor(p(j),:)        
      enddo
      write(*,'(a7,a)')'[---] ','Easiest way is to harcode'
      write(*,'(a7,a)')'[---] ','Look for the edge being falsely detected as bnd'
      write(*,'(a7,a)')'[---] ','Hardocde to remove that edge in bnd'
      write(*,*)
      stop
    else
      do j=1,nEleBnd
        nbnd2=nbnd2+1
        if(nbnd2.gt.nbnd)then
          write(*,'(a7,a)')'[ERR] ','Check bnd detect algortihm loc1'
          write(*,'(a7,a,I10)')'[---] ',' Element = ',i
          write(*,'(a7,a,I10)')'[---] ',' nbnd = ',nbnd
          write(*,'(a7,a,I10)')'[---] ','nbnd2 = ',nbnd2
          stop
        endif
        bnd(nbnd2,:)=eleBnd(j,:)
      enddo      
    endif
  enddo
  if(nbnd2.ne.nbnd)then
    write(*,'(a7,a)')'[ERR] ','Check bnd detect algortihm loc2'
    write(*,'(a7,a,I10)')'[---] ',' nbnd = ',nbnd
    write(*,'(a7,a,I10)')'[---] ','nbnd2 = ',nbnd2
    stop
  endif


  !! Convert to Lon Lat  
  coor(:,1)=(coor(:,1)-sc1cx)*scdx+sc1lon
  coor(:,2)=(coor(:,2)-sc1cy)*scdy+sc1lat
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
  ! write(fl3,'("Number of Nodes")')
  ! write(fl3,'(i15)')np
  ! write(fl3,'("Nodes")')
  do i=1,np
    write(fl3,'(2f20.6)')coor(i,:)
  enddo

  !! Mesh file
  write(fl2,'(3a15)')'Elements','Nodes'
  write(fl2,'(3i15)')nele,np
  write(fl2,'(2a15)')'BndSides','BndTypes'
  write(fl2,'(2i15)')nbnd,nbndty
  write(fl2,'("Nodes")')
  do i=1,np
    write(fl2,'(2f20.6)')coor(i,:)
  enddo  
  write(fl2,'("Elements")')
  do i=1,nele
    write(fl2,'(4i15)')conn(i,:)
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
  
  write(outname,'("Mesh_",I4.4,4I2.2)')timeVal(1:3),timeVal(5:6)
  call out4NXML(outname,np,np,nele,5,0,conn,&
    coor(:,1),coor(:,2),coor(:,1),coor(:,2),coor(:,1),coor(:,2))

  close(fl1)
  close(fl2)
  close(fl3)
  write(*,'(a7,a)')"[MSG] ","Completed successfully"
  201 format(a7,a20,a)
end program quadmesh2DToBsnq

subroutine out4NXML(probname,npoinl,npoint,nelem,code,ts,&
  conn,coorx,coory,u,v,eta,dep)
implicit none
  
  integer(kind=4),intent(in)::npoinl,npoint,nelem
  integer(kind=4),intent(in)::conn(nelem,4),code,ts
  integer(kind=4)::i,j,k

  real(kind=8),intent(in)::coorx(npoint),coory(npoint)
  real(kind=8),intent(in)::u(npoint),v(npoint),eta(npoint)
  real(kind=8),intent(in)::dep(npoint)  

  character(len=200),intent(in)::probname
  character(len=200)::text


  write(text,'(I10)')ts
  text=adjustl(text)
  text=probname(1:len_trim(probname))//"_"//text(1:len_trim(text))//".vtu"
  open(code,file=text(1:len_trim(text)))
  write(code,'(a)')'<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'
  write(code,'(T3,a)')'<UnstructuredGrid>'
  write(code,'(T5,a,i10,a,i10,a)')'<Piece NumberOfPoints="',npoinl,'" NumberOfCells="',nelem,'">'


  write(code,'(T5,a)')'<PointData Scalars="eta" Vectors="vel">'
  
  write(code,'(T7,a)')'<DataArray type="Float64" Name="eta" format="ascii">'
  write(code,*)eta(1:npoinl)
  write(code,'(T7,a)')'</DataArray>'

  write(code,'(T7,a)')'<DataArray type="Float64" Name="depth" format="ascii">'
  write(code,*)dep(1:npoinl)
  write(code,'(T7,a)')'</DataArray>'

  ! write(code,'(T7,a)')'<DataArray type="Float64" Name="porH" format="ascii">'
  ! write(code,*)porH-depth(1:npoinl)
  ! write(code,'(T7,a)')'</DataArray>'

  write(code,'(T7,a)')'<DataArray type="Float64" Name="vel" NumberOfComponents="3" format="ascii">'  
  do i=1,npoinl
    write(code,*)u(i),v(i),0
  enddo
  write(code,'(T7,a)')'</DataArray>'
  
  write(code,'(T5,a)')'</PointData>'  

  write(code,'(T5,a)')'<CellData>'
  write(code,'(T5,a)')'</CellData>'

  write(code,'(T5,a)')'<Points>'
  write(code,'(T7,a)')'<DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">'  
  do i=1,npoinl
    write(code,*)coorx(i),coory(i),0
  enddo
  write(code,'(T7,a)')'</DataArray>'
  write(code,'(T5,a)')'</Points>'

  write(code,'(T5,a)')'<Cells>'
  write(code,'(T7,a)')'<DataArray type="Float64" Name="connectivity" format="ascii">'  
  do i=1,nelem
    write(code,*)conn(i,1:4)-1
  enddo
  write(code,'(T7,a)')'</DataArray>'
  write(code,'(T7,a)')'<DataArray type="Float64" Name="offsets" format="ascii">'  
  do i=1,nelem
    write(code,*)4*i
  enddo
  write(code,'(T7,a)')'</DataArray>'
  write(code,'(T7,a)')'<DataArray type="UInt8" Name="types" format="ascii">'  
  do i=1,nelem
    write(code,*)9
  enddo
  write(code,'(T7,a)')'</DataArray>'
  write(code,'(T5,a)')'</Cells>'

  write(code,'(T5,a)')'</Piece>'
  write(code,'(T3,a)')'</UnstructuredGrid>'
  write(code,'(a)')'</VTKFile>'
  close(code)

end subroutine out4NXML