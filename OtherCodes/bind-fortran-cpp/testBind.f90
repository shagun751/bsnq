program testBind
use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE, C_CHAR, C_NULL_PTR
  
  interface
    function GetObject( n ) result( obj ) BIND(C, name="GetObject")      
      use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE, C_CHAR, C_NULL_PTR

      integer(C_INT), value, intent(in) :: n
      type(C_PTR) :: obj
    end function GetObject

    subroutine display( obj, n ) BIND(C)
      use, intrinsic :: ISO_C_BINDING, only : C_INT, C_PTR, C_DOUBLE, C_CHAR, C_NULL_PTR

      type(C_PTR), value, intent(in) :: obj
      integer(C_INT), value, intent(in) :: n
    end subroutine display
  end interface

  type(C_PTR), save :: obj = C_NULL_PTR
  
  obj = GetObject( 10 )

  call display( obj, 20 )


end program testBind