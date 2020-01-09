!!------------------------basicVars------------------------!!
module basicVars  
implicit none
  
  integer,parameter::C_K1=4,C_K2=8
  integer(kind=C_K1)::i,i1,i2,j,j1,j2,k,k1,k2
  integer(kind=C_K1),parameter::mafi(10)=(/ (i,i=11,20)  /)
  integer(kind=C_K1)::iel,n1,n2,n3,n4,n5,n6
  integer(kind=C_K1)::nq(6),nl(3)
  integer(kind=C_K1)::tmpi1,tmpi2,tmpi3,tmpi4,tmpi5
  integer(kind=8)::sysClk(5),sysRaInt
  
  
  real(kind=C_K2)::tmpr1,tmpr2,tmpr3,tmpr4,tmpr5 
  real(kind=C_K2)::sysRate,sysTime(5)        

  !!----------------------Constants-----------------------!!
  !! Number of Gauss-Quadrature points
  integer(kind=C_K1),parameter::nGP=16
  !! Neighbour nodes and elements limit
  integer(kind=C_K1),parameter::maxNePoi=30
  integer(kind=C_K1),parameter::maxNeEle=10
  !! Gravity
  real(kind=C_K2),parameter::grav=9.81d0          
  !! PI
  real(kind=C_K2),parameter::pi=atan(1d0)*4d0       
  !! Water density
  real(kind=C_K2),parameter::rhoW=1000d0  
  !! Bsnq constant
  real(kind=C_K2),parameter::BsqC=1d0/15d0
  !!--------------------End Constants---------------------!!


  !!  mafi defintions
  !!  mafi(1)     Mesh File
  !!  mafi(2)     Paraview output
  !!  mafi(3)     Volume output
  !!  mafi(4)     <Unknown>
  !!  mafi(5)     Input file
  !!  mafi(6)     Porosity file
  !!  mafi(7)     Wave probes files


end module basicVars
!!----------------------End basicVars----------------------!!

!!------------------------shapeFnc-------------------------!!
module shapeFnc
use basicVars
implicit none
  
  real(kind=C_K2),protected::gpE(nGP),gpN(nGP),gpW(nGP)  
  real(kind=C_K2),protected::sh3F(3,nGP),sh6F(6,nGP)
  real(kind=C_K2),protected::sh3FE(3,nGP),sh6FE(6,nGP)
  real(kind=C_K2),protected::sh3FN(3,nGP),sh6FN(6,nGP)  

  type jacbType
    real(kind=C_K2)::D(nGP),D11(nGP),D12(nGP)
    real(kind=C_K2)::D21(nGP),D22(nGP)    
  end type jacbType

  contains
    !!-----------initGaussPoi-----------!!
    subroutine initGaussPoi    
    implicit none
      if(nGP.eq.16)then
        gpE(1)=1d0/3d0
        gpE(2)=0.45929258829272d0
        gpE(3)=0.45929258829272d0
        gpE(4)=0.08141482341455d0
        gpE(5)=0.17056930775176d0
        gpE(6)=0.17056930775176d0
        gpE(7)=0.65886138449648d0        
        gpE(8)=0.05054722831703d0
        gpE(9)=0.05054722831703d0
        gpE(10)=0.89890554336594d0
        gpE(11)=0.26311282963464d0
        gpE(12)=0.72849239295540d0
        gpE(13)=0.00839477740996d0
        gpE(14)=0.72849239295540d0
        gpE(15)=0.26311282963464d0
        gpE(16)=0.00839477740996d0

        gpN(1)=1d0/3d0
        gpN(2)=0.45929258829272d0
        gpN(3)=0.08141482341455d0        
        gpN(4)=0.45929258829272d0
        gpN(5)=0.17056930775176d0
        gpN(6)=0.65886138449648d0
        gpN(7)=0.17056930775176d0
        gpN(8)=0.05054722831703d0
        gpN(9)=0.89890554336594d0
        gpN(10)=0.05054722831703d0
        gpN(11)=0.72849239295540d0
        gpN(12)=0.00839477740996d0
        gpN(13)=0.26311282963464d0
        gpN(14)=0.26311282963464d0
        gpN(15)=0.00839477740996d0
        gpN(16)=0.72849239295540d0

        gpW(1)=0.14431560767779d0
        gpW(2)=0.09509163426728d0
        gpW(3)=0.09509163426728d0
        gpW(4)=0.09509163426728d0
        gpW(5)=0.10321737053472d0
        gpW(6)=0.10321737053472d0
        gpW(7)=0.10321737053472d0
        gpW(8)=0.03245849762320d0
        gpW(9)=0.03245849762320d0
        gpW(10)=0.03245849762320d0
        gpW(11)=0.02723031417443d0
        gpW(12)=0.02723031417443d0
        gpW(13)=0.02723031417443d0
        gpW(14)=0.02723031417443d0
        gpW(15)=0.02723031417443d0
        gpW(16)=0.02723031417443d0
        gpW=gpW/2d0
      endif
    end subroutine initGaussPoi
    !!---------End initGaussPoi---------!!

    !!-----------initShapeFnc-----------!!
    subroutine initShapeFnc    
    implicit none
      
      real(kind=C_K2)::pe,pn
      
      do i=1,nGP
        pe=gpE(i)
        pn=gpN(i)

        sh6F(1,i)=1d0-3d0*pe+2d0*pe*pe-3d0*pn &
                  +2d0*pn*pn+4d0*pe*pn
        sh6F(2,i)=-pe+2d0*pe*pe
        sh6F(3,i)=-pn+2d0*pn*pn
        sh6F(4,i)=4d0*pe-4d0*pe*pe-4d0*pe*pn
        sh6F(5,i)=4d0*pe*pn;
        sh6F(6,i)=4d0*pn-4d0*pe*pn-4d0*pn*pn

        sh6FE(1,i)=-3d0+4d0*pe+4d0*pn
        sh6FE(2,i)=-1d0+4d0*pe
        sh6FE(3,i)=0d0        
        sh6FE(4,i)=4d0-8d0*pe-4d0*pn
        sh6FE(5,i)=4d0*pn
        sh6FE(6,i)=-4d0*pn

        sh6FN(1,i)=-3d0+4d0*pe+4d0*pn
        sh6FN(2,i)=0d0
        sh6FN(3,i)=-1d0+4d0*pn;
        sh6FN(4,i)=-4d0*pe
        sh6FN(5,i)=4d0*pe
        sh6FN(6,i)=4d0-4d0*pe-8d0*pn

        sh3F(1,i)=1d0-pe-pn
        sh3F(2,i)=pe
        sh3F(3,i)=pn

        sh3FE(1,i)=-1d0
        sh3FE(2,i)=1d0
        sh3FE(3,i)=0d0

        sh3FN(1,i)=-1d0
        sh3FN(2,i)=0d0
        sh3FN(3,i)=1d0
      enddo

    end subroutine initShapeFnc
    !!---------End initShapeFnc---------!!

end module
!!----------------------End shapeFnc-----------------------!!
