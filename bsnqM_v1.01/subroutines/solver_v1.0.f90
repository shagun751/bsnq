subroutine solveSys(np,nnz,iv,jv,gA,gB,gX,errLim,maxiter,&
  iter,resnorm,ier)
use bsnqGlobVars

  interface
    subroutine paralution_fortran_solve_coo( n, m, nnz, solver, & 
      mformat, preconditioner, pformat, &
      rows, cols, rval, rhs, atol, rtol, div, maxiter, basis, &
      p, q, x, iter, resnorm, ierr ) BIND(C)

      use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE, C_CHAR

      integer(kind=C_INT), value, intent(in)  :: n, m, nnz, maxiter, basis, p, q
      real(kind=C_DOUBLE), value, intent(in)  :: atol, rtol, div
      integer(kind=C_INT),        intent(out) :: iter, ierr
      real(kind=C_DOUBLE),        intent(out) :: resnorm
      type(C_PTR),         value, intent(in)  :: rows, cols, rval, rhs
      type(C_PTR),         value              :: x
      character(kind=C_CHAR)                  :: solver, mformat, preconditioner, pformat

    end subroutine paralution_fortran_solve_coo

    subroutine paralution_fortran_solve_csr( n, m, nnz, solver, &
      mformat, preconditioner, pformat, &
      ivCSR, jvCSR, rval, rhs, atol, rtol, div, maxiter, basis, &
      p, q, x, iter, resnorm, ierr ) BIND(C)

      use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE, C_CHAR

      integer(kind=C_INT), value, intent(in)  :: n, m, nnz, maxiter, basis, p, q
      real(kind=C_DOUBLE), value, intent(in)  :: atol, rtol, div
      integer(kind=C_INT),        intent(out) :: iter, ierr
      real(kind=C_DOUBLE),        intent(out) :: resnorm
      type(C_PTR),         value, intent(in)  :: ivCSR, jvCSR, rval, rhs
      type(C_PTR),         value              :: x
      character(kind=C_CHAR)                  :: solver, mformat, preconditioner, pformat

    end subroutine paralution_fortran_solve_csr
  end interface
  
  integer(kind=C_K1),intent(in)::np,nnz,maxIter
  integer(kind=C_K1),intent(in),target::iv(np+1),jv(nnz)
  integer(kind=C_K1),intent(out)::iter,ier
  real(kind=C_K2),intent(in)::errLim  
  real(kind=C_K2),intent(in),target::gA(nnz),gB(np)  
  real(kind=C_K2),intent(out)::resnorm
  real(kind=C_K2),intent(out),target::gX(np)

  call paralution_fortran_solve_csr(np,np,nnz,&
    'GMRES' // C_NULL_CHAR, &
    'CSR' // C_NULL_CHAR, &
    'None' // C_NULL_CHAR, &
    'CSR' // C_NULL_CHAR, &
    C_LOC(iv), C_LOC(jv), &
    C_LOC(gA), C_LOC(gB), &
    errLim, 1e-15_C_DOUBLE, 1e+8_C_DOUBLE, maxIter, &
    30, 0, 1, C_LOC(gX), iter, resnorm, ier)

end subroutine solveSys