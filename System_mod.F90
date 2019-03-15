!! Author:  Christian Schwermann
!! E-mail:  c.schwermann@wwu.de
!! Date:    15/03/2019
!! Project: PEREMB 
!! File:    System_mod.F90 
!! Copyright: © 2019 Christian Schwermann, ALL RIGHTS RESERVED
!!
!! Definition of the GLOBAL total system, which can contain any number of molecules.
module System_mod
   use Output_mod ! Error and warning
   use Types_mod  ! molecule type
   implicit none

   ! Number of molecules in the system, system (array of molecules)
   public :: nsys, system
   ! System initialization
   public :: Sys_init
   ! values hould be saved, so that subroutines can modify the system globally
   save

   !TODO Maybe add a nice subroutine to extend a system 
   !(by creating a larger system and copying everything)
   !! Number of subsystems, the only thing you have to know beforehand
   integer                       :: nsys
   !! THE system, just an array of molecules, but GLOBAL!
   type(molecule_t), allocatable :: system(:)


contains


   !! Basic initialization of the system, not specifying any molecules
   subroutine Sys_init( numsys )
      !! Input: number of subsystems
      integer, intent(in) :: numsys
      !!

      nsys = numsys
      allocate( system(nsys) )

   end subroutine Sys_init

end module System_mod
