!! Author:  Christian Schwermann
!! E-mail:  c.schwermann@wwu.de
!! Date:    15/03/2019
!! Project: PEREMB 
!! File:    System_mod.F90 
!! Copyright: © 2019 Christian Schwermann, ALL RIGHTS RESERVED
!!
!!***********************************************************************
!!
!! Definition of the GLOBAL total system, which can contain any number of molecules.
!!
!!***********************************************************************
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
   subroutine Init( numsys )
      !! Input: number of subsystems
      integer, intent(in) :: numsys
      !!

      nsys = numsys
      allocate( system(nsys) )

   end subroutine Init


   !! Initialize a subsystem, optionally with all important parameters.
   subroutine Sys_init( subsystem, density, vnuc, gridweights, gridpositions, active, spin )
      !! Input: index of the subsystem
      integer, intent(in) :: subsystem
      !! Input: density, 1D array of length ngpt (or 2*ngpt if spin polarized; order is a, b) 
      real(kind=DP), intent(in)       :: density(:)
      !! Input: nuclear potential, of length ngpt
      real(kind=DP), optional, intent(in)       :: vnuc(:)
      !! Input: grid weights, length ngpt
      real(kind=DP), optional, intent(in)       :: gridweights(:)
      !! Input: grid positions, length (3,ngpt)
      real(kind=DP), optional, intent(in)       :: gridpositions(:,:)
      !! Input: whether the system is active
      logical, optional, intent(in) :: active
      !! Input: whether the system is spin polarized
      logical, optional, intent(in) :: spin
      !! Internal :: number of grid points
      integer                         :: ngpt
      !!

      if( .not. Allocated( system ) ) &
         & call Error( "Sys_init: system not initialized!" )

      ngpt = Size( density )
      if( Present( spin ) ) then
         if( spin ) then
            ngpt = ngpt / 2
         end if
      end if

      system(subsystem) = molecule_t(ngpt = ngpt)

      ! set spin first, as density depends on it
      if( Present( spin ) ) &
         & system(subsystem)%spinpol = spin

      call Mol_set_density( system(subsystem), density )

      if( Present( vnuc ) ) &
         & call Mol_set_vnuc( system(subsystem), vnuc )

      if( Present( gridweights ) ) then
         if( Present( gridpositions ) ) then
            call Mol_init_grid( system(subsystem), weights = gridweights, positions = gridpositions )
         else
            call Mol_init_grid( system(subsystem), weights = gridweights )
         end if
      end if

      if( Present( active ) ) &
         & system(subsystem)%active = active

   end subroutine Sys_init


end module System_mod
