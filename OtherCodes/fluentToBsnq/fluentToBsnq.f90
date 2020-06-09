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
  integer(kind=4)::tmpia(3)
  real(kind=8),allocatable::coor(:,:)
  character(len=200)::probname,text,outname
  character(len=100)::tempstr(3)  

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

  call date_and_time(VALUES=timeVal)
  write(outname,'("Mesh_",I4.4,4I2.2,".plt")')timeVal(1:3),timeVal(5:6)
  !write(outname,'("Mesh_",I4.4,4I2.2,".plt")')timeVal(1:3),timeVal(5),2
  write(*,201)"[INF] ",'Output Name: ',trim(outname)
  fl2=11
  open(fl2,file=trim(outname))

!!--------------------------Reading---------------------------!!
  read(fl1,*)text
  read(fl1,*)nele,np,nbnd,nbndty,nedg
  write(*,*)nele,np,nbnd,nbndty,nedg

  !! Mesh file
  write(fl2,'(3a15)')'Elements','Nodes','Edges'
  write(fl2,'(3i15)')nele,np,nedg
  write(fl2,'(2a15)')'BndSides','BndTypes'
  write(fl2,'(2i15)')nbnd,nbndty
  
  read(fl1,*)text
  write(fl2,'("Nodes")')
  allocate(coor(np,2))
  do i=1,np
    read(fl1,*)coor(i,1:2)
    write(fl2,'(2f20.6)')coor(i,:)    
  enddo
  write(*,*)'[MSG] Done node read'

  read(fl1,*)text
  write(fl2,'("Elements")')
  allocate(conn(nele,3))
  do i=1,nele
    read(fl1,*)tempstr(3),tempstr(1),tempstr(2)
    do j=1,3
      call hex2dec(tempstr(j),len_trim(tempstr(j)),conn(i,j),&
        l2)
      if(l2.eq.1)then
        write(*,*)'[ERR] Error in hex2dec'
      endif
    enddo
    write(fl2,'(3i15)')conn(i,:)    
  enddo
  write(*,*)'[MSG] Done element read'

  k=0;
  do l=1,nbndty
    read(fl1,*)text
    read(fl1,*)i,j
    write(fl2,*)trim(text)
    write(fl2,'(2i15)')i,j
    do i=k+1,k+j
      read(fl1,*)(tempstr(j),j=1,3)
      do j=1,3
        call hex2dec(tempstr(j),len_trim(tempstr(j)),tmpia(j),&
          l2)
        if(l2.eq.1)then
          write(*,*)'[ERR] Error in hex2dec'
        endif
      enddo     
      write(fl2,'(3i15)')tmpia(:)
    enddo    
  enddo  
  write(*,*) "[MSG] Boundaries read done" 

  write(fl2,*)'Depth'

  
!!------------------------End Reading-------------------------!!

  close(fl1)
  close(fl2)
  ! close(fl3)
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

subroutine hex2dec(text,varlen,result,err)
implicit none
  integer(kind=4),intent(out)::result,err
  integer(kind=4),intent(in)::varlen
  character(len=1),intent(in)::text(varlen)
  !local
  integer(kind=4)::i,j
  character(len=16)::hexstring

  err=0
  hexstring="0123456789abcdef"
  result=0;
  do i=varlen,1,-1
    j=SCAN(hexstring,text(i))-1    
    if(j.lt.0) then
      err=1
    endif
    result=result+(j*16**(varlen-i))
  enddo
end subroutine hex2dec