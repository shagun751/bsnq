!!---------------------------vertVelMod----------------------------!!
module vertVelMod
use bsnqGlobVars
implicit none

  type, public :: vertVelDerv      
    real(kind=C_K2)::u,ux,uxx,uxxx
    real(kind=C_K2)::uh,uhx,uhxx,uhxxx
    real(kind=C_K2)::hx
    ! d2(U)/dxdt and d2(Uh)/dxdt
    real(kind=C_K2)::uxt,uhxt      
    ! d(U)/dx and d(Uh)/dx at t(n-1), t(n-2)   
    real(kind=C_K2)::ux_tn(2),uhx_tn(2)    
  contains
  !   procedure ::  init => initVertVelDerv
    procedure ::  initByInterp
  end type vertVelDerv

contains

! !!-----------------------initVertVelDerv-----------------------!!
!   subroutine initVertVelDerv(b,np)
!   use bsnqGlobVars
!   implicit none

!     class(vertVelDerv),intent(inout)::b
!     integer(kind=C_K1),intent(in)::np

!     allocate(b%u(np), b%ux(np), b%uxx(np), b%uxxx(np))
!     allocate(b%uh(np), b%uhx(np), b%uhxx(np), b%uhxxx(np))
!     allocate(b%hx(np))
!     allocate(b%uxt(np), b%uhxt(np))
!     allocate(b%uxtn(np,2), b%uhxtn(np,2))
!     b%ux=0d0
!     b%uhx=0d0    
!     b%uxtn=0d0
!     b%uhxtn=0d0    
!     b%np=np    

!   end subroutine initVertVelDerv
! !!---------------------End initVertVelDerv---------------------!!



!!------------------------initByInterp-------------------------!!
  subroutine initByInterp(b,nNei,nei,wei,i)
  implicit none

    class(vertVelDerv),intent(out)::b
    integer(kind=C_K1),intent(in)::nNei
    real(kind=C_K2),intent(in)::wei(nNei)
    type(vertVelDerv),intent(in)::nei(nNei)

    integer(kind=C_K1),intent(out)::i

    b%u=0d0
    b%ux=0d0
    b%uxx=0d0
    b%uxxx=0d0
    b%uh=0d0
    b%uhx=0d0
    b%uhxx=0d0
    b%uhxxx=0d0
    b%hx=0d0
    b%uxt=0d0
    b%uhxt=0d0    

    do i=1,nNei
      b%u = b%u + wei(i) * nei(i)%u
      b%ux =  b%ux + wei(i) * nei(i)%ux
      b%uxx = b%uxx +  wei(i) * nei(i)%uxx
      b%uxxx =  b%uxxx + wei(i) * nei(i)%uxxx
      b%uh =  b%uh + wei(i) * nei(i)%uh
      b%uhx = b%uhx +  wei(i) * nei(i)%uhx
      b%uhxx =  b%uhxx + wei(i) * nei(i)%uhxx
      b%uhxxx = b%uhxxx +  wei(i) * nei(i)%uhxxx
      b%hx =  b%hx + wei(i) * nei(i)%hx
      b%uxt = b%uxt +  wei(i) * nei(i)%uxt
      b%uhxt =  b%uhxt + wei(i) * nei(i)%uhxt      
    enddo

  end subroutine initByInterp
!!----------------------End initByInterp-----------------------!!



!!-------------------------getVertVel--------------------------!!
  subroutine vertVelExp(z,h,eta,hx,u,ux,uxx,uxxx,uhx,uhxx,uhxxx,&
    uxt,uhxt,uc,wc,pc)
  implicit none

    real(kind=C_K2),intent(in)::z,h,eta,hx,u,ux,uxx,uxxx
    real(kind=C_K2),intent(in)::uhx,uhxx,uhxxx,uxt,uhxt
    real(kind=C_K2),intent(out)::uc,wc,pc

    uc = u - (0.5d0*h*uhxx - h*h*uxx/6d0) - (z*uhxx + 0.5d0*z*z*uxx)
    wc = -uhx - z*ux + z/2d0*( hx*uhxx + h*uhxxx ) &
      - z/6d0*( 2d0*h*hx*uxx + h*h*uxxx ) &
      + z*z/2d0*uhxxx + z*z*z/6d0*uxxx
    pc = rhoW * ( grav*( -z + eta ) + z*(uhxt + 0.5d0*z*uxt) )    

  end subroutine vertVelExp
!!-----------------------End getVertVel------------------------!!

end module vertVelMod
!!---------------------------End shipMod---------------------------!!