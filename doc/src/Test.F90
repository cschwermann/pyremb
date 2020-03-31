!! Author:  Christian Schwermann
!! E-mail:  c.schwermann@wwu.de
!! Date:    15/03/2019
!! Project: PEREMB 
!! File:    Test.F90 
!! Copyright: © 2019 Christian Schwermann, ALL RIGHTS RESERVED
!!
!!***********************************************************************
!!
!! Test of most of the functions and subroutines, just to see if they work / are correct
!!
!!***********************************************************************
program Test
   use Precision_mod ! IEEE Precision
   use Types_mod     ! Molecule, grid and ions type
   use System_mod    ! global system
   use Utils_mod     ! Utilities 
   use XCpot_mod     ! DFT 
   implicit none
   !! number of grid points, number of atoms, loop index
   integer :: ngpt, natoms, i
   !! density, nuclear potential, ion charges, ion positions, cell, another density
   real(kind=DP)  :: dens(1:3), vnuc(1:3), q(1:4), pos(1:3, 1:4), cell(1:3, 1:3), dens2(1:2)
   !! three molecules
   type(molecule_t) :: mol1, mol2, mol3
   !! grid
   type(grid_t) :: grid1
   !! dft energy and effective potential
   real(kind=DP) :: ener, pot(1:3)
   !! grid positions for two grids, reference and result function for interpolation test
   real(kind=DP) :: positions1(1:3, 1:4), positions2(1:3, 1:5), func(1:4), res_func(1:5)
   !! 

   write(*,*) "Test Program"

   ! initialization
   ngpt = 3
   dens = (/ 1.0_DP, 2.0_DP, 3.0_DP /)
   vnuc = (/ 4.0_DP, 5.0_DP, 6.0_DP /)
   natoms = 4
   q = (/ 1.0_DP, 2.0_DP, 3.0_DP, 4.0_DP /)
   pos=0.0_DP
   do i = 1, 4
      pos(1, i) = i * 0.1_DP
      pos(2, i) = i * 0.2_DP
      pos(3, i) = i * 0.3_DP
   end do

   cell=0.0_DP
   do i = 1, 3
      cell(i, i) = 5.0_DP
   end do

   ! Test molecule type
   mol1 = molecule_t( ngpt = ngpt, density = dens )
   mol2 = mol1

   ! Test molecule initialization
   call Mol_set_vnuc( mol1, vnuc )
   call Mol_init_ions( mol1, natoms, q, pos )
   call Mol_init_grid( mol1, cell = cell )

   ! Output: "Has"-functions
   write(*,*) Mol_has_vnuc( mol1 ), Mol_has_ions( mol1 ), Mol_has_grid( mol1 )
   write(*,*) mol1%vnuc
   write(*,*) Ions_has_charges( mol1%ions )
   write(*,*) mol1%grid%cell
   write(*,*) Mol_has_vnuc( mol2 ), Mol_has_ions( mol2 ), Mol_has_grid( mol2 )

   ! Test molecule Settings: ions, grid
   call Mol_set_ions( mol2, mol1%ions )
   call Mol_set_grid( mol2, mol1%grid )

   write(*,*) Mol_has_vnuc( mol2 ), Mol_has_ions( mol2 ), Mol_has_grid( mol2 )

   ! Test molecule Settings: spin, density, gradient
   call Mol_set_spin( mol2, .true. )
   call Mol_set_density( mol2, (/ dens, dens /) )
   call Mol_set_gradient( mol2, (/ dens, dens, dens /) )

   ! Different initialization, also works
   mol3 = molecule_t( ngpt = ngpt, spinpol = .true., density = (/ dens, dens /), gradient = (/ dens, dens, dens /) )

   mol3%vnuc = (/ 1.0_DP, 2.0_DP, 3.0_DP, 4.0_DP /)

   ! Output: "Has"-functions
   write(*,*) Mol_has_density( mol3 ), Mol_has_gradient( mol3 ), Mol_has_grid( mol3 ), Mol_has_vnuc( mol3 )

   ! Test system nitialization
   call Sys_init( 3 )

   system(1) = mol1
   system(2) = mol2
   system(3) = mol3
   
   system(2)%active = .true.

   write(*,*) "System:", nsys, system(1)%active, system(2)%active, system(3)%active

   ! Test integration
   write(*,*) "Density integral:", Integrate( dens )

   grid1=grid_T( & 
      & ngpt = 3, &
      & weights = (/ 0.1_DP, 0.2_DP, 0.3_DP /) )

   write(*,*) "Weighted density integral:", Integrate( dens, grid1 )

   ! Test DFT procedures for libxc
   call Xc_energy( mol1, "LDA_XC_TETER93", ener )
   write(*,*) "LDA Energy:", ener

   call Xc_potential( mol1, "LDA_XC_TETER93", pot )
   write(*,*) "LDA Potential:", pot

   call Xc_energy( mol3, "GGA_XC_PBE1W", ener )
   write(*,*) "PBE Energy:", ener
   
   positions1(:, :) = 0.0_DP
   positions1(1, :) = (/ 1.0_DP, 3.0_DP, 5.0_DP, 7.0_DP /)
   positions2(:, :) = 0.0_DP
   positions2(1, :) = (/ 0.0_DP, 2.0_DP, 4.0_DP, 6.0_DP, 8.0_DP /)
   ! simple x**2 function as test
   func(:) = positions1(1, :) ** 2

   ! Test Shepard interpolation
   call Shepard_interpolate( positions1, func, positions2, res_func )

   ! Print reference
   write(*,*) "Interpolation: reference"
   write(*,*) "Integral:", Integrate( func ) * 2.0_DP
   write(*,'(A5,A5)') "X", "Y"

   do i = 1, 4
      write(*,'(F5.1,F5.1)') positions1(1, i), func(i)
   end do

   ! Print result
   write(*,*) "Interpolation: result"
   write(*,*) "Integral:", Integrate( res_func ) * 2.0_DP
   write(*,'(A5,A5)') "X", "Y"

   do i = 1, 5
      write(*,'(F5.1,F5.1)') positions2(1, i), res_func(i)
   end do

   ! Also see what happens if we try to get our reference from the result
   call Shepard_interpolate( positions2, res_func, positions1, func )

   ! Print reference obtained from interpolating twice
   write(*,*) "Interpolation: reference backwards"
   write(*,*) "Integral:", Integrate( func ) * 2.0_DP
   write(*,'(A5,A5)') "X", "Y"

   do i = 1, 4
      write(*,'(F5.1,F5.1)') positions1(1, i), func(i)
   end do

   ! Also check exponential / sinh weights
   call Shepard_interpolate( positions2, res_func, positions1, func, texp = .true. )

   ! Print results
   write(*,*) "Interpolation: reference backwards, sinh weights"
   write(*,*) "Integral:", Integrate( func ) * 2.0_DP
   write(*,'(A5,A5)') "X", "Y"

   do i = 1, 4
      write(*,'(F5.1,F5.1)') positions1(1, i), func(i)
   end do

   ! Test already known positions
   positions1(1, :) = (/ 5.0_DP, 3.0_DP, 5.0_DP, 7.0_DP /)
   positions2(1, :) = (/ 2.0_DP, 3.0_DP, 5.0_DP, 6.0_DP, 8.0_DP /)

   call Shepard_interpolate( positions1, func, positions2, res_func )

   ! Print results
   write(*,*) "Interpolation: result"
   write(*,*) "Integral:", Integrate( res_func ) * 2.0_DP
   write(*,'(A5,A5)') "X", "Y"

   do i = 1, 5
      write(*,'(F5.1,F5.1)') positions2(1, i), res_func(i)
   end do

   !TODO test Mol_interpolate

   ! Test error termination
   dens2 = (/ 1.0_DP, 2.0_DP /)
   call Mol_set_density( mol2, dens2 )

end program Test
