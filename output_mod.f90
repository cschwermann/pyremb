!Output Module
!contains error and warning subroutine
module Output_mod
   use Precision_mod !consistent precisions
   implicit none

   public :: Error, Warning

contains

   !Throw error and exit with -1
   subroutine Error( string, int1, int2, real1, real2 )
      character(len=*), intent(in) :: string
      integer, optional, intent(in) :: int1, int2
      real(kind=DP), optional, intent(in) :: real1, real2

      write(*,*) "ERROR! ", string
      if( Present( int1 ) ) write(*,*) "(1)=", int1
      if( Present( int2 ) ) write(*,*) "(2)=", int2
      if( Present( real1 ) ) write(*,*) "(1)=", real1
      if( Present( real2 ) ) write(*,*) "(2)=", real2
      stop
   end subroutine Error

   !Throw warning, don't exit!
   subroutine Warning( string, int1, int2, real1, real2 )
      character(len=*), intent(in) :: string
      integer, optional, intent(in) :: int1, int2
      real(kind=DP), optional, intent(in) :: real1, real2

      write(*,*) "WARNING! ", string
      if( Present( int1 ) ) write(*,*) "(1)=", int1
      if( Present( int2 ) ) write(*,*) "(2)=", int2
      if( Present( real1 ) ) write(*,*) "(1)=", real1
      if( Present( real2 ) ) write(*,*) "(2)=", real2
   end subroutine Warning

end module Output_mod
