program fnpt2bsnq
use fnptModule
implicit none
  
  integer::i
  double precision:: t,eta,p,q,t0

  call readFNPT
  open(9,file='output.dat')

  call depIntVelFNPT(t0,eta,p)    

  q=0d0
  do while(.not.fnptEOF)
    
    call depIntVelFNPT(t,eta,p)    
    write(9,'(4F20.8)')t-t0,eta,p,q
    call updateFNPT    
  enddo
  close(9)

end program fnpt2bsnq
