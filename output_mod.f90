!Output Module
!contains error subroutine
module output_mod
   implicit none

 contains

   !Throw error and exit with -1
   subroutine error(string,int1,int2,real1,real2)
      character(*),intent(in) :: string
      integer,optional,intent(in) :: int1,int2
      real*8,optional,intent(in) :: real1,real2

      write(*,*) "ERROR! ",string
      if(present(int1)) write(*,*) "(1)=",int1
      if(present(int2)) write(*,*) "(2)=",int2
      if(present(real1)) write(*,*) "(1)=",real1
      if(present(real2)) write(*,*) "(2)=",real2
      stop
   end subroutine

   !Throw warning, don't exit!
   subroutine warning(string,int1,int2,real1,real2)
      character(*),intent(in) :: string
      integer,optional,intent(in) :: int1,int2
      real*8,optional,intent(in) :: real1,real2

      write(*,*) "WARNING! ",string
      if(present(int1)) write(*,*) "(1)=",int1
      if(present(int2)) write(*,*) "(2)=",int2
      if(present(real1)) write(*,*) "(1)=",real1
      if(present(real2)) write(*,*) "(2)=",real2
   end subroutine

end module
