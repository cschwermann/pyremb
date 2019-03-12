!Defines the GLOBAL total system, which can contain any number of molecules.
module System_mod
   use Output_mod
   use Types_mod
   implicit none

   public :: nsys, system
   public :: Sys_init
   save

   !Number of systems, the only thing you have to know beforehand
   !TODO Maybe add a nice subroutine to extend a system 
   !(by creating a larger system and copying everything)
   integer                      :: nsys
   !THE system, just an array of molecules, but GLOBAL!
   type(molecule_t), allocatable, dimension(:) :: system

contains

   !basic initialization of the system, not specifying any molecules
   subroutine Sys_init( numsys )
      integer, intent(in) :: numsys

      nsys = numsys
      allocate( system(nsys) )

   end subroutine Sys_init

end module System_mod
