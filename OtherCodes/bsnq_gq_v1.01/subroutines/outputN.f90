subroutine output(npoinl,npoint,nelem,conn,coord,probname,mafi,&
  rtime,p,q,eta,w,depth,pdt,qdt,etaMin,etaMax,porH)
implicit none
  
  integer(kind=4),intent(in)::npoinl,npoint,nelem
  integer(kind=4),intent(in)::conn(nelem,6),mafi(10)
  integer(kind=4)::i,j,k,code

  real(kind=8),intent(in)::coord(npoint,2),p(npoint),w(npoinl)
  real(kind=8),intent(in)::q(npoint),eta(npoinl),depth(npoint)
  real(kind=8),intent(in)::pdt(npoint),qdt(npoint),rtime
  real(kind=8),intent(in)::etaMin(npoinl),etaMax(npoinl)
  real(kind=8),intent(in)::porH(npoinl)

  character(len=100),intent(in)::probname
  character(len=100)::text

  code=mafi(2)

  write(text,'(I15)')int(rtime*1000000)
  text=adjustl(text)
  text="Output/"//probname(1:len_trim(probname))//"_"//text(1:len_trim(text))//".vtk"
  open(code,file=text(1:len_trim(text)))
  write(code,'(a26)')'# vtk DataFile Version 3.0'
  write(code,'(a35)')'Instananeous Variables Output'
  write(code,'(a5)')'ASCII'
  write(code,'(a25)')'DATASET UNSTRUCTURED_GRID'
  write(code,'(a6,1x,i10,1x,a6)')'POINTS',npoinl,'double'
  do i=1,npoinl
    write(code,*)coord(i,1),coord(i,2),0d0
  enddo

  !nelem*(nep+1)
  write(code,'(a5,1x,i10,1x,i10)')'CELLS',nelem,3*nelem+nelem 
  do i=1,nelem
    write(code,'(500(i10,1x))')3,(conn(i,j)-1,j=1,3)
  enddo
  write(code,'(a10,1x,i10)')'CELL_TYPES',nelem
  do i=1,nelem
   write(code,'(i1)')5
  enddo

  write(code,'(a10,1x,i10)')'POINT_DATA',npoinl
  write(code,*)'VECTORS ', "velocity ", "float "
  do i=1,npoinl
    write(code,*)p(i),q(i),0.0d0
  enddo

  ! write(code,*)'VECTORS ', "accel ", "float "
  ! do i=1,npoinl
  !   write(code,*)pdt(i),qdt(i),0.0d0
  ! enddo

  write(code,*)'SCALARS ', "eta ", "float ",1
  write(code,'(a20)')'LOOKUP_TABLE default'
  do i=1,npoinl
    write(code,*)eta(i)
  enddo

  write(code,*)'SCALARS ', "waveHeight ", "float ",1
  write(code,'(a20)')'LOOKUP_TABLE default'
  do i=1,npoinl
    write(code,*)etaMax(i)-etaMin(i)
  enddo

  write(code,*)'SCALARS ', "depth ", "float ",1
  write(code,'(a20)')'LOOKUP_TABLE default'
  do i=1,npoinl
    write(code,*)-depth(i)
  enddo

  write(code,*)'SCALARS ', "porH ", "float ",1
  write(code,'(a20)')'LOOKUP_TABLE default'
  do i=1,npoinl
    write(code,*)porH(i)-depth(i)
  enddo

  close(code)

end subroutine output

subroutine outputPoiXML(npoinl,npoint,nelem,conn,coord,probname,mafi,&
  rtime,p,q,eta,w,depth,pdt,qdt,etaMin,etaMax,porH)
implicit none
  
  integer(kind=4),intent(in)::npoinl,npoint,nelem
  integer(kind=4),intent(in)::conn(nelem,6),mafi(10)
  integer(kind=4),allocatable::tmp1(:),tmp2(:)
  integer(kind=4)::i,j,k,code

  real(kind=8),intent(in)::coord(npoint,2),p(npoint),w(npoinl)
  real(kind=8),intent(in)::q(npoint),eta(npoinl),depth(npoint)
  real(kind=8),intent(in)::pdt(npoint),qdt(npoint),rtime
  real(kind=8),intent(in)::etaMin(npoinl),etaMax(npoinl)
  real(kind=8),intent(in)::porH(npoinl)

  character(len=100),intent(in)::probname
  character(len=100)::text

  allocate(tmp1(npoinl),tmp2(npoinl))
  do i=1,npoinl
    tmp1(i)=i
  enddo
  tmp2=1  

  code=mafi(2)

  write(text,'(I15)')int(rtime*1000000)
  text=adjustl(text)
  text="Output/"//probname(1:len_trim(probname))//"_"//text(1:len_trim(text))//".vtu"
  open(code,file=text(1:len_trim(text)))
  write(code,'(a)')'<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'
  write(code,'(T3,a)')'<UnstructuredGrid>'
  write(code,'(T5,a,i10,a,i10,a)')'<Piece NumberOfPoints="',npoinl,'" NumberOfCells="',npoinl,'">'

  write(code,'(T5,a)')'<PointData Scalars="eta" Vectors="vel">'

  write(code,'(T7,a)')'<DataArray type="Float64" Name="eta" format="ascii">'
  write(code,*)eta
  write(code,'(T7,a)')'</DataArray>'

  write(code,'(T7,a)')'<DataArray type="Float64" Name="waveH" format="ascii">'
  write(code,*)etaMax-etaMin
  write(code,'(T7,a)')'</DataArray>'

  write(code,'(T7,a)')'<DataArray type="Float64" Name="depth" format="ascii">'
  write(code,*)-depth(1:npoinl)
  write(code,'(T7,a)')'</DataArray>'

  write(code,'(T7,a)')'<DataArray type="Float64" Name="porH" format="ascii">'
  write(code,*)porH-depth(1:npoinl)
  write(code,'(T7,a)')'</DataArray>'

  write(code,'(T7,a)')'<DataArray type="Float64" Name="vel" NumberOfComponents="3" format="ascii">'  
  do i=1,npoinl
    write(code,*)p(i),q(i),0
  enddo
  write(code,'(T7,a)')'</DataArray>'

  write(code,'(T5,a)')'</PointData>'  

  write(code,'(T5,a)')'<CellData>'
  write(code,'(T5,a)')'</CellData>'

  write(code,'(T5,a)')'<Points>'
  write(code,'(T7,a)')'<DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">'  
  ! ! write(code,'(T5,a)',advance='no')'  '
  ! write(code,*)transpose(coord(1:10,:))
  do i=1,npoinl
    write(code,*)coord(i,:),0
  enddo
  write(code,'(T7,a)')'</DataArray>'
  write(code,'(T5,a)')'</Points>'

  write(code,'(T5,a)')'<Cells>'
  write(code,'(T7,a)')'<DataArray type="Float64" Name="connectivity" format="ascii">'  
  write(code,*)tmp1-1
  write(code,'(T7,a)')'</DataArray>'
  write(code,'(T7,a)')'<DataArray type="Float64" Name="offsets" format="ascii">'  
  write(code,*)tmp1
  write(code,'(T7,a)')'</DataArray>'
  write(code,'(T7,a)')'<DataArray type="UInt8" Name="types" format="ascii">'  
  write(code,*)tmp2
  write(code,'(T7,a)')'</DataArray>'
  write(code,'(T5,a)')'</Cells>'

  write(code,'(T5,a)')'</Piece>'
  write(code,'(T3,a)')'</UnstructuredGrid>'
  write(code,'(a)')'</VTKFile>'
  close(code)
  deallocate(tmp1,tmp2)

end subroutine outputPoiXML

subroutine outputXML(npoinl,npoint,nelem,conn,coord,probname,&
  mafi,rtime,p,q,eta,w,depth)
implicit none
  
  integer(kind=4),intent(in)::npoinl,npoint,nelem
  integer(kind=4),intent(in)::conn(nelem,6),mafi(10)
  integer(kind=4),allocatable::tmp1(:),tmp2(:)
  integer(kind=4)::i,j,k,code

  real(kind=8),intent(in)::coord(npoint,2),p(npoint),w(npoinl)
  real(kind=8),intent(in)::q(npoint),eta(npoinl),depth(npoint)
  real(kind=8),intent(in)::rtime
  ! real(kind=8),intent(in)::etaMin(npoinl),etaMax(npoinl)  

  character(len=100),intent(in)::probname
  character(len=100)::text

  code=mafi(2)

  write(text,'(I15)')int(rtime*1000000)
  text=adjustl(text)
  text="Output/"//probname(1:len_trim(probname))//"_"//text(1:len_trim(text))//".vtu"
  open(code,file=text(1:len_trim(text)))
  write(code,'(a)')'<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'
  write(code,'(T3,a)')'<UnstructuredGrid>'
  write(code,'(T5,a,i10,a,i10,a)')'<Piece NumberOfPoints="',npoinl,'" NumberOfCells="',nelem,'">'

  write(code,'(T5,a)')'<PointData Scalars="eta" Vectors="vel">'
  
  write(code,'(T7,a)')'<DataArray type="Float64" Name="eta" format="ascii">'
  write(code,*)eta
  write(code,'(T7,a)')'</DataArray>'

  ! write(code,'(T7,a)')'<DataArray type="Float64" Name="waveH" format="ascii">'
  ! write(code,*)etaMax-etaMin
  ! write(code,'(T7,a)')'</DataArray>'

  write(code,'(T7,a)')'<DataArray type="Float64" Name="depth" format="ascii">'
  write(code,*)-depth(1:npoinl)
  write(code,'(T7,a)')'</DataArray>'

  ! write(code,'(T7,a)')'<DataArray type="Float64" Name="porH" format="ascii">'
  ! write(code,*)porH-depth(1:npoinl)
  ! write(code,'(T7,a)')'</DataArray>'

  write(code,'(T7,a)')'<DataArray type="Float64" Name="vel" NumberOfComponents="3" format="ascii">'  
  do i=1,npoinl
    write(code,*)p(i),q(i),0
  enddo
  write(code,'(T7,a)')'</DataArray>'
  
  write(code,'(T5,a)')'</PointData>'  

  write(code,'(T5,a)')'<CellData>'
  write(code,'(T5,a)')'</CellData>'

  write(code,'(T5,a)')'<Points>'
  write(code,'(T7,a)')'<DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">'  
  do i=1,npoinl
    write(code,*)coord(i,:),0
  enddo
  write(code,'(T7,a)')'</DataArray>'
  write(code,'(T5,a)')'</Points>'

  write(code,'(T5,a)')'<Cells>'
  write(code,'(T7,a)')'<DataArray type="Float64" Name="connectivity" format="ascii">'  
  do i=1,nelem
    write(code,*)conn(i,1:3)-1
  enddo
  write(code,'(T7,a)')'</DataArray>'
  write(code,'(T7,a)')'<DataArray type="Float64" Name="offsets" format="ascii">'  
  do i=1,nelem
    write(code,*)3*i
  enddo
  write(code,'(T7,a)')'</DataArray>'
  write(code,'(T7,a)')'<DataArray type="UInt8" Name="types" format="ascii">'  
  do i=1,nelem
    write(code,*)5
  enddo
  write(code,'(T7,a)')'</DataArray>'
  write(code,'(T5,a)')'</Cells>'

  write(code,'(T5,a)')'</Piece>'
  write(code,'(T3,a)')'</UnstructuredGrid>'
  write(code,'(a)')'</VTKFile>'
  close(code)

end subroutine outputXML
