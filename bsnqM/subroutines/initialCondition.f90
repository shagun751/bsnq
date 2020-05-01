!! Apply initial conditions inside the subroutine initMat() inside 
!! bsnqModule.f90

!!-----------------------------solitIC-----------------------------!!
  subroutine solitIC(npl,npt,cor,er,pr,qr)
  use bsnqGlobVars
  implicit none
    
    integer(kind=C_K1),intent(in)::npl,npt
    real(kind=C_K2),intent(in)::cor(npt,2)
    real(kind=C_K2),intent(out)::er(npl),pr(npt),qr(npt)

    integer(kind=C_K1)::i
    real(kind=C_K2)::tmpr1,tmpr2,tmpr3,tmpr4,tmpr5

    ! pr=0d0
    ! qr=0d0
    ! !er=0.045d0*dexp(-2d0*( (cor(1:npl,1)-18.288d0)**2 ))
    ! do i=1,npl
    !   er(i)=(cor(i,1)-18.288d0)**2
    !   if(er(i).gt.25)then
    !     er(i)=0d0
    !   else
    !     er(i)=0.045*exp(-2d0*er(i))
    !   endif
    ! enddo

    er=0d0
    pr=0d0
    qr=0d0
    tmpr1=0.45d0
    tmpr2=0.045d0
    tmpr3=dsqrt(grav*(tmpr1+tmpr2))
    tmpr4=dsqrt(3*tmpr2/(4*(tmpr1**3)))
    do i=1,npt
      if((cor(i,1).ge.3d0).and.(cor(i,1).le.19d0)) then
        tmpr5=tmpr2/(dcosh(tmpr4*(cor(i,1)-(11d0)))**2)
        pr(i)=tmpr3*tmpr5
        if(i.le.npl) then
          er(i)=tmpr5
        endif
      endif
    enddo  

  end subroutine solitIC
!!---------------------------End solitIC---------------------------!!