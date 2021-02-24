!!-------------------- Version M 1.x.x ------------------------!!
!!  -> Quadratic + Linear
!!  -> Boundary - SemiDirect + Penalty + Gauss Seidel
!!    -> 11 - Inlet - No Absorbing
!!    -> 12 - NoSlip Wall 
!!    -> 13 - Slip Wall - Rectangular wall only
!!    -> 14 - Sponge - BOUSS2D approach generalised input
!!    -> 15 - Outlet  - Not coded
!!  -> Porosity - Emergent structure only (not done yet)
!!    -> Generalised input
!!  -> Solver - Normalised
!!  -> Generalised input - v3
!!  -> GMRES
!!  -> Dirichlet BndCond for eta, p, q
!!  -> Paralution CSR
!!	-> XML output
!!-----------------------------------------------------------------!!
!! Time-Stepping : RK4 

program boussinesqQuad
use bsnqGlobVars
use bsnqModule
implicit none

!!--------------------------Declarations---------------------------!!
  type(bsnqCase)::bq  
  character(len=C_KSTR)::bqtxt

  real(kind=C_K2),allocatable::u(:),v(:),eta(:)
  real(kind=C_K2),allocatable::sedNelX(:),sedNelY(:)
  real(kind=C_K2),allocatable::sedVanX(:),sedVanY(:)
!!------------------------End Declarations-------------------------!!
  
  call getarg(1,bq%probname)  
  do while(len_trim(bq%probname).lt.1)
    write(*,'(A)',advance='no')"Enter Problem Name:"
    read(*,*)bq%probname
  enddo
  write(*,*)"Problem Name: "//trim(bq%probname)

  bqtxt=trim(bq%probname)//'.rout'
  open(9,file=trim(bqtxt))    

  call system_clock(bq%sysC(1))
  call bq%meshRead
  call bq%femInit
  call bq%setRun 
  call bq%setMFree 
  call bq%initMat  
  call bq%statMatrices    

  allocate( u(bq%npt), v(bq%npt), eta(bq%npt) )
  allocate( sedNelX(bq%npt), sedNelY(bq%npt) )
  allocate( sedVanX(bq%npt), sedVanY(bq%npt) )

  call bq%caseOutputs
  
  do while(abs(bq%tOb(0)%rtm-bq%endTime).gt.bq%dt/2d0)
  
    call bq%preInstructs

    call bq%timeStepRK4
  
    call bq%postInstructs   

    u = bq%tOb(0)%p / bq%tOb(0)%tD
    v = bq%tOb(0)%q / bq%tOb(0)%tD
    eta = bq%tOb(0)%tD - bq%dep 

    call sedimentTrans(bq, bq%npt, bq%dep, bq%cor(:,1), bq%cor(:,2), &
      u, v, eta, sedNelX, sedNelY, sedVanX, sedVanY)

    bq%sedNelX = sedNelX
    bq%sedNelY = sedNelY
    bq%sedVanX = sedVanX
    bq%sedVanY = sedVanY

    call bq%caseOutputs

  enddo

  call system_clock(bq%sysC(2))
  close(bq%wpEle(-1))
  write(9,*)"[MSG] boussinesqQuad End"
  write(9,'(" [TIM] ",F15.4)')1d0*(bq%sysC(2)-bq%sysC(1))/bq%sysRate
  close(9)  
end program boussinesqQuad



subroutine sedimentTrans(bq, npt, dep, corx, cory, u, v, eta, &
  sedNelX, sedNelY, sedVanX, sedVanY)
use bsnqGlobVars
use bsnqModule
implicit none
  
  type(bsnqCase),intent(in):: bq
  integer(kind=C_K1),intent(in):: npt
  real(kind=C_K2),intent(in):: dep(npt), u(npt), v(npt), eta(npt)
  real(kind=C_K2),intent(in):: corx(npt), cory(npt)
  real(kind=C_K2),intent(out):: sedNelX(npt), sedNelY(npt)
  real(kind=C_K2),intent(out):: sedVanX(npt), sedVanY(npt)

  !! Note: 
  !! Do not allocate the above variables again.
  !! They have already been allocated

  integer(kind=C_K1):: i
  !real(kind=C_K2):: tmpr
  real(kind=C_K2) :: d50,d90
  real(kind=C_K2) :: rho,nu,dstar,z0,rhos,z0d,fca
  real(kind=C_K2) :: vonkc,gra,s
  real(kind=C_K2), allocatable :: z0_x(:),z0_y(:),depth(:)
  real(kind=C_K2), allocatable ::thetas_x(:),thetas_y(:)
  real(kind=C_K2), allocatable ::tauc_x(:),tauc_y(:)
  real(kind=C_K2),allocatable :: drag_coef(:)
  real(kind=C_K2),allocatable :: Ts_x(:),Ts_y(:),lamdas(:)
  real(kind=C_K2),allocatable :: deltas_x(:),deltas_y(:)
  real(kind=C_K2),allocatable :: z0s_x(:),z0s_y(:),z0t(:)
  real(kind=C_K2) :: thetacr
  real(kind=C_K2) :: lamdar,deltar,ar,taucr,z0r
  real(kind=C_K2) :: beta,value,phi,phiw,phic,phised,bedslope
  real(kind=C_K2), allocatable ::Q_x(:),Q_y(:),bl_x(:),bl_y(:)
  real(kind=C_K2),allocatable :: drag_coef_x(:),drag_coef_y(:)
  real(kind=C_K2),allocatable :: ucr(:),b_x(:),b_y(:),a_x(:),a_y(:) 
  real(kind=C_K2),allocatable :: phy_x(:),phy_y(:)
  real(kind=C_K2),allocatable :: gfk_x(:),hfk_x(:),dfk_y(:),efk_y(:),nyk_x(:),myk_y(:)

  real(kind=C_K2),allocatable :: sedNelX_dX(:), sedNelX_dY(:)
 
!*************************************************************************
!
!              CONSTANTS
!
!*************************************************************************
!
!  ^ y
!  |
!  |
!  |
!  |_  phi (taken clockwise
!  | -->    from y axis)
!  |   /  
!  |  /           
!  | /
!  |_______________> x
!


!   constants & inputs
      
      !pi=4.*datan(1d0)
    !  twopi = 2.*pi
      vonkc=0.4
      gra=9.81
      rhos=2650.                     
      phised = 0.5585                !angle of repose of sediment(32deg)
      rho=1027
      nu=1.36e-06
      d50=4e-04
      d90=4e-04
      phiw=0.0
      phic=0.0
      bedslope=0.0 
      s = rhos/rho             
      z0d=d50/12.    !skin friction roughness
	   			
      allocate(thetas_x(npt),thetas_y(npt),tauc_x(npt),tauc_y(npt))
      allocate(Q_x(npt),Q_y(npt),bl_x(npt),bl_y(npt))
      allocate(z0_x(npt),z0_y(npt),depth(npt))
      allocate(drag_coef(npt))
      allocate(Ts_x(npt),Ts_y(npt),lamdas(npt))	 
      allocate(deltas_x(npt),deltas_y(npt))	 
      allocate(z0s_x(npt),z0s_y(npt),z0t(npt))
      allocate(drag_coef_x(npt),drag_coef_y(npt))
      allocate(ucr(npt),b_x(npt),b_y(npt),a_x(npt),a_y(npt))
      allocate(phy_x(npt),phy_y(npt))
      allocate(gfk_x(npt),hfk_x(npt),dfk_y(npt),efk_y(npt),nyk_x(npt),myk_y(npt))

      allocate( sedNelX_dX(npt), sedNelX_dY(npt) )
      
      phi = phiw-phic


  !***************************************************************************
!
! Compute change in the water level
!
!***************************************************************************
	
	do i=1,npt
         depth(i)=dep(i)+eta(i)
    end do
      
!***************************************************************************
      
!***************************************************************************      
!      
! Calculate Bed Shear Stresses for currents
!    
!
!***************************************************************************

!	uses z0=z0d (i.e. grain size roughness)

     z0=z0d
      
     do i=1,npt
	   
	   drag_coef(i)=(vonkc/(1.+dlog(z0/depth(i))))**2
       
	   ! in x direction
	   
	   tauc_x(i)=rho*drag_coef(i)*u(i)**2
       
	   ! in y direction
	   
	   tauc_y(i)=rho*drag_coef(i)*v(i)**2

       end do
  
      do i=1,npt
	  
	  thetas_x(i) = tauc_x(i)/(gra*(rhos-rho)*d50)
	  thetas_y(i) = tauc_y(i)/(gra*(rhos-rho)*d50)

     end do

      dstar =real (((gra*(s-1)/(nu **2))**(1.0/3.0))*d50)
      thetacr = (0.3/(1+(1.2*dstar)))+(0.055*(1-exp(-0.02*dstar)))

	  ! adjust shield's no for bedslope effects
      beta=ATAN(bedslope)
      value=(COS(beta))**2*(TAN(phised))**2-(SIN(phi))**2*(SIN(beta))**2
      thetacr=thetacr*(COS(phi)*SIN(beta)+SQRT(value))/TAN(phised)
      
! JAS 09/11/04 
! calculate ripple roughness length 

	   do i=1,npt
		 
	 if (thetas_x(i) .GE. thetacr) THEN
       
               taucr = gra*(rhos-rho)*d50*thetacr
               lamdar = 1000*d50
               deltar = lamdar/7
               lamdas(i) = 7.3*depth(i)
               Ts_x(i) = (tauc_x(i)-taucr)/taucr
               deltas_x(i) = 0.11*depth(i)*(d50/depth(i))**(0.3)*(1-exp(-0.5*Ts_x(i)))*(25-Ts_x(i))
               ar = 1.0
               z0r = (ar*(deltar)**2)/lamdar
               z0s_x(i) = (ar*(deltas_x(i))**2)/lamdas(i)
           
           if (thetas_x(i) .LE. 0.8) THEN
               
               z0_x(i) = z0r+z0s_x(i)+z0d
           
             endif
             
           if ( thetas_x(i) .GT. 0.8) then
             
             z0t(i) = (5*tauc_x(i))/(30*gra*(rhos-rho))
             z0_x(i) = z0r+z0s_x(i)+z0d+z0t(i)
             
           endif
           
           else 
           z0_x(i) = z0d
      
      endif 

	  if (thetas_y(i) .GE. thetacr) THEN
       
       taucr = gra*(rhos-rho)*d50*thetacr
               lamdar = 1000*d50
               deltar = lamdar/7
               lamdas(i) = 7.3*depth(i)
               Ts_y(i) = (tauc_y(i)-taucr)/taucr
               deltas_y(i) = 0.11*depth(i)*(d50/depth(i))**(0.3)*(1-exp(-0.5*Ts_y(i)))*(25-Ts_y(i))
               ar = 1.0
               z0r = (ar*(deltar)**2)/lamdar
               z0s_y(i) = (ar*(deltas_y(i))**2)/lamdas(i)
           
           if (thetas_y(i) .LE. 0.8) THEN
               
               z0_y(i) = z0r+z0s_y(i)+z0d
           
             endif
             
           if ( thetas_y(i) .GT. 0.8) then
             
             z0t(i) = (5.0*tauc_y(i))/(30*gra*(rhos-rho))
             z0_y(i) = z0r+z0s_y(i)+z0d+z0t(i)
             
           endif

         else 
           z0_y(i) = z0d
      
      endif 
      
      end do
	  

! calculate total stresses using new roughness length
! will be the same as skin stresses if z0=z0d

      do i=1,npt

	   drag_coef_x(i)=(vonkc/(1.+dlog(z0_x(i)/depth(i))))**2
	   drag_coef_y(i)=(vonkc/(1.+dlog(z0_y(i)/depth(i))))**2
           tauc_x(i)=rho*drag_coef_x(i)*u(i)**2
           tauc_y(i)=rho*drag_coef_y(i)*v(i)**2

	   end do
         !  open(18,file='z0_x.txt', Status = 'Unknown')  

     ! write(18,*) z0_x
         
    
       do i=1,npt
	  
	  thetas_x(i) = tauc_x(i)/(gra*(rhos-rho)*d50)
	  thetas_y(i) = tauc_y(i)/(gra*(rhos-rho)*d50)

     end do
      
      
!**************************************************************************
!
! Calculate bed load (kg/s)
!
!**************************************************************************       
  
    
 ! By Neilson 
      
     do i=1,npt

  ! in x-direction
      
	   if (thetas_x(i) .LT. thetacr) then
       
         !  print*, 'no sediment transport by bedload in x direction'
            sedNelx(i) = 0.0
       
      endif
      
      if (thetas_x(i) .GE. thetacr) THEN
     
         phy_x(i) = 12*(thetas_x(i)**0.5)*(thetas_x(i)-thetacr)                                                     
         bl_x(i) = ((gra*(s-1)*(d50**3.0))**0.5)*phy_x(i) 
         sedNelX(i) = bl_x(i)
		! bl_x(i) = bl_x(i)*rhos                                                
         
		 !PRINT *,'BEDLOAD(kg/m/s) by Neilson IN X-DIRECTION IS = ',bl_x
    
      endif

      if (u(i) .LT. 0) then

      sedNelX(i) = -bl_x(i)

      endif

   ! in y-direction

      if (thetas_y(i) .LT. thetacr) then
       
         !  print*, 'no sediment transport by bedload in y direction'
            sedNelY(i) = 0.0      
     
       endif
      
      if (thetas_y(i) .GE. thetacr) THEN
     
         phy_y(i) = 12*(thetas_y(i)**0.5)*(thetas_y(i)-thetacr)                                                     
         bl_y(i) = ((gra*(s-1)*(d50**3.0))**0.5)*phy_y(i) 
         sedNelY(i) = bl_y(i) 
		! bl_y(i) = bl_y(i)*rhos                                                
         
		! PRINT *,'BEDLOAD(kg/m/s) by Neilson IN Y-DIRECTION IS = ',bl_y
    
      endif

      if (v(i) .LT. 0) then

      sedNelY(i) = -bl_y(i)

      endif

	  end do
	  

    call bq%getAllPoiGrad(sedNelX, sedNelX_dX, 1)
    call bq%getAllPoiGrad(sedNelX, sedNelX_dY, 2)
  !*********************************************************************************************
     sedVanX = 0.0
     sedVanY = 0.0

end subroutine sedimentTrans
!**************************************************************************************************************************
