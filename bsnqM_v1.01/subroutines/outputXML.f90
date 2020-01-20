!!----------------------------outputXML----------------------------!!
  subroutine outputXML(b)
  implicit none

    class(bsnqCase),intent(in)::b 

    write(bqtxt,'(I15)')int(b%tOb(0)%rtm*1000)
    bqtxt=adjustl(bqtxt)
    bqtxt="Output/"//trim(b%probname)//"_"//trim(bqtxt)//".vtu"
    open(newunit=mf,file=trim(bqtxt))

    write(mf,'(a)')'<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'
    write(mf,'(T3,a)')'<UnstructuredGrid>'
    write(mf,'(T5,a,i10,a,i10,a)')'<Piece NumberOfPoints="',b%npl,'" NumberOfCells="',b%nele,'">'

    ! PointData
    write(mf,'(T5,a)')'<PointData Scalars="eta" Vectors="vel">'

    write(mf,'(T7,a)')'<DataArray type="Float64" Name="eta" format="ascii">'
    write(mf,'(F20.6)')b%tOb(0)%e
    write(mf,'(T7,a)')'</DataArray>'

    write(mf,'(T7,a)')'<DataArray type="Float64" Name="waveH" format="ascii">'
    write(mf,'(F20.6)')b%etaMax-b%etaMin
    write(mf,'(T7,a)')'</DataArray>'

    write(mf,'(T7,a)')'<DataArray type="Float64" Name="press" format="ascii">'
    write(mf,'(F15.6)')b%presr(1:b%npl)
    write(mf,'(T7,a)')'</DataArray>'

    ! write(mf,'(T7,a)')'<DataArray type="Float64" Name="absC" format="ascii">'
    ! write(mf,*)b%absC(1:b%npl)
    ! write(mf,'(T7,a)')'</DataArray>'

    write(mf,'(T7,a)')'<DataArray type="Float64" Name="depth" format="ascii">'
    write(mf,'(F15.6)')-b%dep(1:b%npl)
    write(mf,'(T7,a)')'</DataArray>'

    ! write(mf,'(T7,a)')'<DataArray type="Float64" Name="porH" format="ascii">'
    ! write(mf,*)porH-depth(1:npoinl)
    ! write(mf,'(T7,a)')'</DataArray>'

    write(mf,'(T7,a)')'<DataArray type="Float64" Name="vel" NumberOfComponents="3" format="ascii">'  
    do i=1,b%npl
      write(mf,'(2F20.6,F5.1)')b%tOb(0)%p(i), b%tOb(0)%q(i), 0d0
    enddo
    write(mf,'(T7,a)')'</DataArray>'

    write(mf,'(T5,a)')'</PointData>'


    ! CellData
    write(mf,'(T5,a)')'<CellData>'
    write(mf,'(T5,a)')'</CellData>'


    ! Mesh
    write(mf,'(T5,a)')'<Points>'
    write(mf,'(T7,a)')'<DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">'  
    do i=1,b%npl
      write(mf,'(2F15.4,F5.1)')b%cor(i,:),0d0
    enddo
    write(mf,'(T7,a)')'</DataArray>'
    write(mf,'(T5,a)')'</Points>'

    write(mf,'(T5,a)')'<Cells>'
    write(mf,'(T7,a)')'<DataArray type="Int32" Name="connectivity" format="ascii">'  
    do i=1,b%nele
      write(mf,'(3I12)')b%conn(i,1:3)-1
    enddo
    write(mf,'(T7,a)')'</DataArray>'
    write(mf,'(T7,a)')'<DataArray type="Int32" Name="offsets" format="ascii">'  
    do i=1,b%nele
      write(mf,'(I12)')3*i
    enddo
    write(mf,'(T7,a)')'</DataArray>'
    write(mf,'(T7,a)')'<DataArray type="UInt8" Name="types" format="ascii">'  
    do i=1,b%nele
      write(mf,'(I12)')5
    enddo
    write(mf,'(T7,a)')'</DataArray>'
    write(mf,'(T5,a)')'</Cells>'

    write(mf,'(T5,a)')'</Piece>'
    write(mf,'(T3,a)')'</UnstructuredGrid>'
    write(mf,'(a)')'</VTKFile>'
    close(mf)

    write(9,'(" [MSG] OutXML at ",F15.4)')b%tOb(0)%rtm

  end subroutine outputXML
!!--------------------------End outputXML--------------------------!!