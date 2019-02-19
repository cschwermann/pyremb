!Defines the GLOBAL total system, which can contain any number of molecules.
module system_mod
   use output_mod
   use types_mod
   implicit none
   public
   save

   !Number of systems, the only thing you have to know beforehand
   !TODO Maybe add a nice subroutine to extend a system 
   !(by creating a larger system and copying everything)
   integer                      :: nsys
   !THE system, just an array of molecules, but GLOBAL!
   type(T_Molecule),allocatable,dimension(:) :: system

 contains

   !basic initialization of the system, not specifying any molecules
   subroutine Sys_init(numsys)
      integer,intent(in) :: numsys
      integer            :: i

      nsys=numsys
      allocate(system(nsys))

   end subroutine


end module
