subroutine resume_input(npoinl,npoint,rtime,p,q,eta,text)
implicit none

  integer(kind=4),intent(in)::npoinl,npoint
  integer(kind=4)::ifile,i,j,k

  real(kind=8),intent(out)::p(npoint),q(npoint),eta(npoinl),rtime

  character(len=100),intent(in)::text
  character(len=200)::text2

  logical::ex

  p=0d0
  q=0d0
  eta=0d0

  ifile=111
  inquire(file=text(1:len_trim(text)),exist=ex)
  if(ex) then
    open(ifile,file=text(1:len_trim(text)))
  else
    write(9,*)"[Error] Input file missing"
    stop
  endif

  read(ifile,*)text2
  read(ifile,*)rtime

  read(ifile,*)text2
  do i=1,npoint
    read(ifile,*)p(i),q(i)
  enddo

  read(ifile,*)text2
  do i=1,npoinl
    read(ifile,*)eta(i)
  enddo

  close(ifile)

  write(9,*)"[Msg] Resume input successfull!"


end subroutine resume_input

subroutine resume_output(npoinl,npoint,probname,mafi,rtime,p,q,eta)
implicit none

  integer(kind=4),intent(in)::npoinl,npoint,mafi(10)
  integer(kind=4)::i,j,k,code

  real(kind=8),intent(in)::p(npoint),q(npoint),eta(npoinl),rtime

  character(len=100),intent(in)::probname
  character(len=100)::text

  code=mafi(4)

  write(text,'(I15)')int(rtime*1000000)
  text=adjustl(text)
  text="Resume/resume_"//probname(1:len_trim(probname))//"_"//text(1:len_trim(text))//".vtk"
  open(code,file=text(1:len_trim(text)))

  write(code,*)"Resume Time"
  write(code,*)rtime
  write(code,*)"velocity vector"
  do i=1,npoint
    write(code,*)p(i),q(i)
  enddo
  write(code,*)"eta scalar"
  do i=1,npoinl
    write(code,*)eta(i)
  enddo

  close(code)

end subroutine resume_output