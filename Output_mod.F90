!! Author:  Christian Schwermann
!! E-mail:  c.schwermann@wwu.de
!! Date:    15/03/2019
!! Project: PEREMB 
!! File:    Output_mod.F90 
!! Copyright: © 2019 Christian Schwermann, ALL RIGHTS RESERVED
!!
!!***********************************************************************
!!
!! Definition of Error and Warning functions used throughout all of PEREMB
!!
!!***********************************************************************
module Output_mod
   use Precision_mod ! IEEE precision
   implicit none

   ! Error and Warning functions
   public :: Error, Warning

contains


   !!********************************************************************
   !! Throw an Error determined by a string, print some numbers and 
   !! exit with non-zero status
   !!********************************************************************
   subroutine Error( string, int1, int2, real1, real2 )
      !! Input: string describing the error
      character(len=*), intent(in) :: string
      !! Input: integers to print out, optional
      integer, optional, intent(in) :: int1, int2
      !! Input: floats to print out, optional
      real(kind=DP), optional, intent(in) :: real1, real2
      !!

      write(*,*) "ERROR! ", string

      if( Present( int1 ) ) write(*,*) "(1)=", int1
      if( Present( int2 ) ) write(*,*) "(2)=", int2
      if( Present( real1 ) ) write(*,*) "(1)=", real1
      if( Present( real2 ) ) write(*,*) "(2)=", real2

      error stop

   end subroutine Error


   !!********************************************************************
   !! Throw a Warnig determined by a string, print some numbers, but
   !! do not exit
   !!********************************************************************
   subroutine Warning( string, int1, int2, real1, real2 )
      !! Input: string describing the error
      character(len=*), intent(in) :: string
      !! Input: integers to print out, optional
      integer, optional, intent(in) :: int1, int2
      !! Input: floats to print out, optional
      real(kind=DP), optional, intent(in) :: real1, real2
      !!

      write(*,*) "WARNING! ", string

      if( Present( int1 ) ) write(*,*) "(1)=", int1
      if( Present( int2 ) ) write(*,*) "(2)=", int2
      if( Present( real1 ) ) write(*,*) "(1)=", real1
      if( Present( real2 ) ) write(*,*) "(2)=", real2

   end subroutine Warning


end module Output_mod
