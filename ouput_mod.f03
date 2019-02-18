!Optput Module
!contains error subroutine
module output_mod
   implicit none

 contains
   !Throw error and exit with -1
   subroutine error(string)
      character(*),intent(in) :: string

      write(*,*) string
      call exit(-1)
   end subroutine

end module
