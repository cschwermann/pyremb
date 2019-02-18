!Types Module
!Contains all used types, including
!Molecule, Grid and Ions
module types
   implicit none
   private

   !Molecule Class, contains points, grid, density, nuclear potential and ions
   type public T_Molecule
      !#Grid Points
      integer               :: ngpt
      !Grid, may contain weights and cell
      type(T_Grid)          :: grid=T_Grid(ngpt=ngpt)
      !Electron Density
      real*8                :: density(ngpt)
      !Nuclear Potential, alternative to ions
      real*8,optional       :: vnuc(ngpt)
      !Ions,contains #ions, charges and positions
      type(T_Ions),optional :: ions=T_Ions()

    contains
      procedure :: set_points    => Mol_set_points
      procedure :: set_grid      => Mol_set_grid
      procedure :: set_density   => Mol_set_density
      procedure :: set_vnuc      => Mol_set_vnuc
      procedure :: set_ions      => Mol_set_ions

      procedure :: get_points    => Mol_get_points
      procedure :: get_grid      => Mol_get_grid
      procedure :: get_density   => Mol_get_density
      procedure :: get_vnuc      => Mol_get_vnuc
      procedure :: get_ions      => Mol_get_ions

      procedure :: has_vnuc      => Mol_has_vnuc
      procedure :: has_ions      => Mol_has_ions
      procedure :: has_grid      => Mol_has_grid
      procedure :: has_density   => Mol_has_density

      procedure :: init          => Mol_init
   end type T_Molecule

   !Grid Class, contains points, weights and cell
   type public T_Grid
      !#Grid Points
      integer               :: ngpt
      !weights,optional
      real*8,optional       :: weights(ngpt)
      !cell vectors,optional
      real*8,optional       :: cell(3,3)

    contains
      procedure :: set_points    => Grid_set_points
      procedure :: set_weights   => Grid_set_weights
      procedure :: set_cell      => Grid_set_cell

      procedure :: get_points    => Grid_get_points
      procedure :: get_weights   => Grid_get_weights
      procedure :: get_cell      => Grid_get_cell

      procedure :: has_weights   => Grid_has_weights
      procedure :: has_cell      => Grid_has_cell
      
      procedure :: init          => Grid_init
   end type T_Grid

   !Ions class, contains number, charges and positions
   type public T_Ions
      !#ions
      integer               :: nions
      !Charges
      real*8                :: charges(nions)
      !Positions,first index is x,y,z
      real*8                :: positions(3,nions)

    contains
      procedure :: set_nions     => Ions_set_nions
      procedure :: set_charges   => Ions_set_charges
      procedure :: set_positions => Ions_set_positions

      procedure :: get_nions     => Ions_get_nions
      procedure :: get_charges   => Ions_get_charges
      procedure :: get_positions => Ions_get_positions

      procedure :: has_charges   => Ions_has_charges

      procedure :: init          => Ions_init
   end type T_Ions

end
