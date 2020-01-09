subroutine waveCalculator(g,pi,errLim,wp,d,wl)
implicit none
  
  integer(kind=4)::i,iterMax

  real(kind=8),intent(inout)::wl
  real(kind=8),intent(in)::wp,d,g,pi,errLim
  real(kind=8)::l0,newl,oldl,x

  iterMax=50000
  l0 = (g/2d0/pi)*wp**2
  oldl = l0  
  do i = 1,iterMax
    newl = l0*dtanh(2d0*pi*d/oldl)
    !if(mod(i,100).eq.0) write(*,*)i,oldl,newl    
    x = abs(newl-oldl)
    if (x.le.errLim) then
      wl = newl
      exit
    else
      oldl = newl
    end if
  end do

  if(i.ge.iterMax) then
    write(*,*)"waveCalculator Error waveL",newl
    wl=-999
  endif

end subroutine wavecalculator