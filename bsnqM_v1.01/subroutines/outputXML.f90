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

    ! write(mf,'(T7,a)')'<DataArray type="Float64" Name="phi" format="ascii">'
    ! do i=1,b%npl
    !   tmpr1=0d0
    !   k=(i-1)* b%ivl(0)
    !   do j=k+1, k+b%ivl(i)
    !     tmpr1=tmpr1 + b%mlsl%phi(j) * b%cor(b%linkl(j),1)
    !   enddo
    !   !if(abs(tmpr1).gt.10) tmpr1=-1d0
    !   write(mf,'(F15.6)')tmpr1
    ! enddo
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

    ! write(mf,'(T7,a)')'<DataArray type="Float64" Name="grad" NumberOfComponents="3" format="ascii">'  
    ! do i=1,b%npl
    !   tmpr1=0d0
    !   tmpr2=0d0
    !   k=(i-1)* b%ivl(0)
    !   do j=k+1, k+b%ivl(i)
    !     tmpr1=tmpr1 + b%mlsl%phiDx(j) * b%cor(b%linkl(j),1)
    !     tmpr2=tmpr2 + b%mlsl%phiDy(j) * b%cor(b%linkl(j),2)
    !   enddo
    !   if(abs(tmpr1).gt.10d0)tmpr1=-1
    !   if(abs(tmpr2).gt.10d0)tmpr2=-1
    !   write(mf,'(2F15.6,F5.1)')tmpr1,tmpr2,0d0
    ! enddo
    ! write(mf,'(T7,a)')'</DataArray>'

    ! write(mf,'(T7,a)')'<DataArray type="Float64" Name="gradEta" NumberOfComponents="3" format="ascii">'  
    ! do i=1,b%npl
    !   tmpr1=0d0
    !   tmpr2=0d0
    !   k=(i-1)* b%ivl(0)
    !   do j=k+1, k+b%ivl(i)
    !     tmpr1=tmpr1 + b%mlsl%phiDx(j) * b%tOb(0)%e(b%linkl(j))
    !     tmpr2=tmpr2 + b%mlsl%phiDy(j) * b%tOb(0)%e(b%linkl(j))
    !   enddo
    !   if(abs(tmpr1).gt.100d0)tmpr1=-10
    !   if(abs(tmpr2).gt.100d0)tmpr2=-10
    !   write(mf,'(2F15.6,F5.1)')tmpr1,tmpr2,0d0
    ! enddo
    ! write(mf,'(T7,a)')'</DataArray>'

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