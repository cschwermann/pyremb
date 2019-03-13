program Test
   use Precision_mod
   use Types_mod
   use System_mod
   use Utils_mod
   use XCpot_mod
   implicit none
   integer :: ngpt, natoms, i
   real(kind=DP)  :: dens(1:3), vnuc(1:3), q(1:4), pos(1:3, 1:4), cell(1:3, 1:3), dens2(1:2)
   type(molecule_t) :: mol1, mol2, mol3
   type(grid_t) :: grid1
   real(kind=DP) :: ener, pot(1:3)

   write(*,*) "Test Program"

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

   mol1 = molecule_t( ngpt = ngpt, density = dens )
   mol2 = mol1

   call Mol_set_vnuc( mol1, vnuc )
   call Mol_init_ions( mol1, natoms, q, pos )
   call Mol_init_grid( mol1, cell = cell )

   write(*,*) Mol_has_vnuc( mol1 ), Mol_has_ions( mol1 ), Mol_has_grid( mol1 )
   write(*,*) mol1%vnuc
   write(*,*) Ions_has_charges( mol1%ions )
   write(*,*) mol1%grid%cell
   write(*,*) Mol_has_vnuc( mol2 ), Mol_has_ions( mol2 ), Mol_has_grid( mol2 )

   call Mol_set_ions( mol2, mol1%ions )
   call Mol_set_grid( mol2, mol1%grid )

   write(*,*) Mol_has_vnuc( mol2 ), Mol_has_ions( mol2 ), Mol_has_grid( mol2 )

   call Mol_set_spin( mol2, .true. )
   call Mol_set_density( mol2, (/ dens, dens /) )
   call Mol_set_gradient( mol2, (/ dens, dens, dens /) )

   mol3 = molecule_t( ngpt = ngpt, spinpol = .true., density = (/ dens, dens /), gradient = (/ dens, dens, dens /) )

   mol3%vnuc = (/ 1.0_DP, 2.0_DP, 3.0_DP, 4.0_DP /)

   write(*,*) Mol_has_density( mol3 ), Mol_has_gradient( mol3 ), Mol_has_grid( mol3 ), Mol_has_vnuc( mol3 )

   call Sys_init( 3 )

   system(1) = mol1
   system(2) = mol2
   system(3) = mol3
   
   system(2)%active = .true.

   write(*,*) "System:", nsys, system(1)%active, system(2)%active, system(3)%active

   write(*,*) "Density integral:", Integrate( dens )

   grid1=grid_T( & 
      & ngpt = 3, &
      & weights = (/ 0.1_DP, 0.2_DP, 0.3_DP /) )

   write(*,*) "Weighted density integral:", Integrate( dens, grid1 )

   call Xc_energy( mol1, "LDA_XC_TETER93", ener )
   write(*,*) "LDA Energy:", ener

   call Xc_potential( mol1, "LDA_XC_TETER93", pot )
   write(*,*) "LDA Potential:", pot

   call Xc_energy( mol3, "GGA_XC_PBE1W", ener )
   write(*,*) "PBE Energy:", ener
   


   dens2 = (/ 1.0_DP, 2.0_DP /)
   call Mol_set_density( mol2, dens2 )

end program Test
