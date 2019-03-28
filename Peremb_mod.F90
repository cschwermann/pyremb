!! Author:  Christian Schwermann
!! E-mail:  c.schwermann@wwu.de
!! Date:    14/03/2019
!! Project: PEREMB 
!! File:    Peremb_mod.F90 
!! Copyright: © 2019 Christian Schwermann, ALL RIGHTS RESERVED
!!
!!***********************************************************************
!!
!! Main subroutines for periodic subystem DFT calculations
!!
!!***********************************************************************
module Peremb_mod
   use Precision_mod ! IEEE precision
   use Output_mod    ! Error and warning
   use utils_mod     ! Some maths
   use System_mod    ! System and subsystems 
   implicit none

   ! Embedding potential for predefined system
   public :: Embedding_potential
   ! Embedding energy for predefined system
   public :: Embedding_energy

contains


   !!********************************************************************
   !! Calculation of embedding potential for the global system.
   !! 
   !! Only needs eXchange, Correlation and Kinetic functionals as input
   !! and returns the effective embedding potential
   !!********************************************************************
   subroutine Embedding_potential( x_func, c_func, k_func, potential )
      !! Input: exchange functional
      character(*), intent(in) :: x_func
      !! Input: correlation functional
      character(*), intent(in) :: c_func
      !! Input: kinetic functional
      character(*), intent(in) :: k_func
      !! Output: embedding potential
      real(kind=DP), intent(out) :: potential(:)
      !! Internal: temporary molecule, total system
      type(molecule_t) :: molecule, total
      !! Internal: temporary grid
      type(grid_t) :: grid
      !! Internal: #points for active system, index of active subsystem, loop index
      integer :: ngpt, iactive, i
      !! Internal: temporary potential
      real(kind=DP), allocatable :: temppot(:)
      !! 

      do i = 1, nsys
         if( system(i)%active ) then
            ngpt = system(i)%ngpt
            grid = system(i)%grid
            total = system(i) 
            iactive = i
         end if
      end do

      potential(1:ngpt) = 0.0_DP
      allocate( temppot(1:ngpt) )

      do i = 1, nsys
         molecule = system(i)
         if( .not. molecule%active ) then
            call Mol_interpolate( molecule, grid )
            if( Mol_has_vnuc( molecule ) ) then
               potential(:) = potential(:) + molecule%vnuc(:)
            else
               !TODO calculate vnuc from ions
               call Error( "Embedding_potential: Need vnuc!" )
            end if
         end if
         call Mol_set_density( total, total%density(:) + molecule%density(:) )
      end do

      ! Total system DFT energy
      call Xc_potential( total, x_func, temppot )
      potential(:) = potential(:) + temppot(:)
      call Xc_potential( total, c_func, temppot )
      potential(:) = potential(:) + temppot(:)
      call Xc_potential( total, k_func, temppot )
      potential(:) = potential(:) + temppot(:)

      ! active system DFT energy
      call Xc_potential( system(iactive), x_func, temppot )
      potential(:) = potential(:) - temppot(:)
      call Xc_potential( system(iactive), c_func, temppot )
      potential(:) = potential(:) - temppot(:)
      call Xc_potential( system(iactive), k_func, temppot )
      potential(:) = potential(:) - temppot(:)

      deallocate( temppot )

   end subroutine Embedding_potential


   !!********************************************************************
   !! Calculation of embedding energy for the global system.
   !! 
   !! Only needs eXchange, Correlation and Kinetic functionals as input
   !! and returns the effective embedding energy on the active grid
   !! @TODO calculate energy on a common supergrid @ENDTODO
   !!********************************************************************
   subroutine Embedding_energy( x_func, c_func, k_func, energy )
      !! Input: exchange functional
      character(*), intent(in) :: x_func
      !! Input: correlation functional
      character(*), intent(in) :: c_func
      !! Input: kinetic functional
      character(*), intent(in) :: k_func
      !! Output: embedding potential
      real(kind=DP), intent(out) :: energy
      !! Internal: temporary molecule
      type(molecule_t) :: molecule, total
      !! Internal: temporary grid
      type(grid_t) :: grid
      !! Internal: #points for active system, index of active system, loop index
      integer :: ngpt, iactive, i
      !! Internal: temporary energy
      real(kind=DP) :: tempener
      !! 

      do i = 1, nsys
         if( system(i)%active ) then
            ngpt = system(i)%ngpt
            grid = system(i)%grid
            total = system(i) 
            iactive = i
         end if
      end do

      energy = 0.0_DP
      do  i = 1, nsys
         molecule = system(i)
         if( .not. molecule%active ) then
            call Mol_interpolate( molecule, grid )
            if( Mol_has_vnuc( molecule ) ) then
               call Mol_set_vnuc( total, total%vnuc(:) + molecule%vnuc(:) )
            else
               !TODO calculate vnuc from ions
               call Error( "Embedding_potential: Need vnuc!" )
            end if
         end if
         call Mol_set_density( total, total%density(:) + molecule%density(:) )
      end do

      ! Electrostatic energy
      energy = energy + Integrate( system(iactive)%vnuc(:) * total%density(:), grid )
      energy = energy + Integrate( molecule%vnuc(:) * system(iactive)%density(:), grid )

      ! Total system DFT energy
      call Xc_energy( total, x_func, tempener )
      energy = energy + tempener
      call Xc_energy( total, c_func, tempener )
      energy = energy + tempener
      call Xc_energy( total, k_func, tempener )
      energy = energy + tempener

      ! active system DFT energy
      call Xc_energy( system(iactive), x_func, tempener )
      energy = energy - tempener
      call Xc_energy( system(iactive), c_func, tempener )
      energy = energy - tempener
      call Xc_energy( system(iactive), k_func, tempener )
      energy = energy - tempener

   end subroutine Embedding_energy

end module Peremb_mod
            
