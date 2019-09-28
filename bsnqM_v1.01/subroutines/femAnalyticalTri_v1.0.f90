subroutine fem_N6i_Sc3_dN3jdx(mat,h1,h2,h3,b11,b12)
use bsnqGlobVars
implicit none

  real(kind=C_K2),intent(out)::mat(6,3)
  real(kind=C_K2),intent(in)::h1,h2,h3,b11,b12

  mat(1,1)=(-(1d0/120d0))*(b11 + b12)*(2d0*h1 - h2 - h3)
  mat(1,2)=(1d0/120d0)*b11*(2d0*h1 - h2 - h3)
  mat(1,3)=(1d0/120d0)*b12*(2d0*h1 - h2 - h3)

  mat(2,1)=(1d0/120d0)*(b11 + b12)*(h1 - 2d0*h2 + h3)
  mat(2,2)=(-(1d0/120d0))*b11*(h1 - 2d0*h2 + h3)
  mat(2,3)=(-(1d0/120d0))*b12*(h1 - 2d0*h2 + h3)

  mat(3,1)=(1d0/120d0)*(b11 + b12)*(h1 + h2 - 2d0*h3)
  mat(3,2)=(-(1d0/120d0))*b11*(h1 + h2 - 2d0*h3)
  mat(3,3)=(-(1d0/120d0))*b12*(h1 + h2 - 2d0*h3)

  mat(4,1)=(-(1d0/30d0))*(b11 + b12)*(2d0*h1 + 2d0*h2 + h3)
  mat(4,2)=(1d0/30d0)*b11*(2d0*h1 + 2d0*h2 + h3)
  mat(4,3)=(1d0/30d0)*b12*(2d0*h1 + 2d0*h2 + h3)

  mat(5,1)=(-(1d0/30d0))*(b11 + b12)*(h1 + 2d0*(h2 + h3))
  mat(5,2)=(1d0/30d0)*b11*(h1 + 2d0*(h2 + h3))
  mat(5,3)=(1d0/30d0)*b12*(h1 + 2d0*(h2 + h3))

  mat(6,1)=(-(1d0/30d0))*(b11 + b12)*(2d0*h1 + h2 + 2d0*h3)
  mat(6,2)=(1d0/30d0)*b11*(2d0*h1 + h2 + 2d0*h3)
  mat(6,3)=(1d0/30d0)*b12*(2d0*h1 + h2 + 2d0*h3)

end subroutine fem_N6i_Sc3_dN3jdx