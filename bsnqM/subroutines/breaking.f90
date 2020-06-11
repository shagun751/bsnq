!!--------------------------breakRoller1---------------------------!!
  subroutine breakRoller1(npt, tD, u, v, p, q, tDavgt1, pObf,&
   brkOn, Rxx)
  use bsnqGlobVars
  use meshFreeMod
  implicit none

    integer(kind=C_K1),intent(in)::npt
    real(kind=C_K2),intent(in)::tD(npt), u(npt), v(npt)
    real(kind=C_K2),intent(in)::p(npt), q(npt), tDavgt1(npt)
    type(mfPoiTyp),intent(in)::pObf(npt)
    real(kind=C_K2),intent(out)::Rxx(npt)
    logical(kind=C_LG),intent(inout)::brkOn(npt)    

    integer(kind=C_K1)::i, j, i2, k2
    real(kind=C_K2)::etadt, c, pDx, qDy, tanPhi
    real(kind=C_K2)::tanPhiB, tanPhi0, rMax, r, del

    Rxx = 0d0    

    tanPhiB = dtan(deg2rad*20)
    tanPhi0 = dtan(deg2rad*10)
    rMax = 0.1638d0

    do i = 1, npt
      i2 = pObf(i)%bsnqId
      k2 = pObf(i)%nn
      call calcGrad(k2, pObf(i)%phiDx, p(pObf(i)%neid), &
        pDx, j)
      call calcGrad(k2, pObf(i)%phiDy, q(pObf(i)%neid), &
        qDy, j)
      c = dsqrt(grav*tD(i))
      etadt = - pDx - qDy

      if(etadt .gt. (c*tanPhiB))then 
        brkOn(i)=.true.
      elseif(etadt .lt. (c*tanPhi0))then
        brkOn(i)=.false.
      endif

      if(brkOn(i).eq..false.) cycle

      tanPhi = etadt / c
      r = min(rMax, 0.45d0*tanPhi)
      del = r*(tD(i) - tDavgt1(i))
      Rxx(i) = del * (1d0 - del/tD(i)) * ((c - u(i))**2)
    enddo

  end subroutine breakRoller1
!!------------------------End breakRoller1-------------------------!!