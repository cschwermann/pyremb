!Types Module
!Contains all used types, including
!Molecule, Grid and Ions
module types_mod
   use output_mod
   implicit none
   private

   !Grid Class, contains points, weights and cell
   type, public :: T_Grid
      !#Grid Points
      integer               :: ngpt
      !weights,optional
      real*8,allocatable    :: weights(:)
      !cell vectors,optional
      real*8,allocatable    :: cell(:,:)

    contains
      procedure :: set_weights   => Grid_set_weights
      procedure :: set_cell      => Grid_set_cell

      procedure :: has_weights   => Grid_has_weights
      procedure :: has_cell      => Grid_has_cell
   end type T_Grid

   !Ions class, contains number, charges and positions
   type, public :: T_Ions
      !#ions
      integer               :: nions
      !Charges
      real*8,allocatable    :: charges(:)
      !Positions,first index is x,y,z
      real*8,allocatable    :: positions(:,:)

    contains
      procedure :: set_charges   => Ions_set_charges
      procedure :: set_positions => Ions_set_positions

      procedure :: has_charges   => Ions_has_charges
      procedure :: has_positions => Ions_has_positions

   end type T_Ions

   !Molecule Class, contains points, grid, density, nuclear potential and ions
   type, public :: T_Molecule
      !#Grid Points
      integer                  :: ngpt
      !Grid, may contain weights and cell
      type(T_Grid),pointer     :: grid => null()
      !Electron Density
      real*8,allocatable       :: density(:)
      !Nuclear Potential, alternative to ions
      real*8,allocatable       :: vnuc(:)
      !Ions,contains #ions, charges and positions
      type(T_Ions),pointer     :: ions => null()

    contains
      procedure :: set_grid      => Mol_set_grid
      procedure :: set_vnuc      => Mol_set_vnuc
      procedure :: set_ions      => Mol_set_ions
      procedure :: set_density   => Mol_set_density

      procedure :: init_grid     => Mol_init_grid
      procedure :: init_ions     => Mol_init_ions

      procedure :: has_vnuc      => Mol_has_vnuc
      procedure :: has_ions      => Mol_has_ions
      procedure :: has_grid      => Mol_has_grid
      procedure :: has_density   => Mol_has_density
   end type T_Molecule

 contains

!Molecule functions

   function Mol_has_grid(this) result(hasgrid)
      Class(T_Molecule),intent(in) :: this
      Logical                      :: hasgrid

      hasgrid = associated(this%grid)
   end function

   function Mol_has_vnuc(this) result(hasvnuc)
      Class(T_Molecule),intent(in) :: this
      Logical                      :: hasvnuc

      hasvnuc = allocated(this%vnuc)
   end function

   function Mol_has_ions(this) result(hasions)
      Class(T_Molecule),intent(in) :: this
      Logical                      :: hasions

      hasions = associated(this%ions)
   end function

   function Mol_has_density(this) result(hasdensity)
      Class(T_Molecule),intent(in) :: this
      Logical                      :: hasdensity

      hasdensity = allocated(this%density)
   end function

   subroutine Mol_set_grid(this,grid)
      Class(T_Molecule),intent(inout) :: this
      type(T_Grid),intent(in)      :: grid

      if(grid%ngpt.ne.this%ngpt) &
         & call error("Mol_set_grid: grid%ngpt /= this%ngpt",grid%ngpt,this%ngpt)

      if(.not.this%has_grid()) allocate(this%grid)
      this%grid=grid
   end subroutine

   subroutine Mol_set_vnuc(this,vnuc)
      Class(T_Molecule),intent(inout) :: this
      real*8,intent(in)      :: vnuc(:)

      if(size(vnuc).ne.this%ngpt) &
         & call error("Mol_set_vnuc: size(vnuc) /= this%ngpt",size(vnuc),this%ngpt)

      if(.not.this%has_vnuc()) allocate(this%vnuc(this%ngpt))
      this%vnuc=vnuc(1:this%ngpt)
   end subroutine

   subroutine Mol_set_ions(this,ions)
      Class(T_Molecule),intent(inout) :: this
      type(T_Ions),intent(in)      :: ions

      if(.not.this%has_ions()) allocate(this%ions)
      this%ions=ions
   end subroutine

   subroutine Mol_set_density(this,density)
      Class(T_Molecule),intent(inout) :: this
      real*8,intent(in)            :: density(:)

      if(size(density).ne.this%ngpt) &
         & call error("Mol_set_density: size(density) /= this%ngpt",size(density),this%ngpt)

      if(.not.this%has_density()) allocate(this%density(this%ngpt))
      this%density=density(1:this%ngpt)
   end subroutine

   subroutine Mol_init_grid(this,weights,cell)
      Class(T_Molecule),intent(inout) :: this
      real*8,optional              :: weights(:)
      real*8,optional              :: cell(:,:)

      if(this%has_grid()) &
         & call error("Mol_init_grid: grid already initialized")

      allocate(this%grid)
      this%grid=T_Grid(ngpt=this%ngpt)

      if(present(weights)) then
         call this%grid%set_weights(weights)
      endif

      if(present(cell)) then
         call this%grid%set_cell(cell)
      endif
   end subroutine

   subroutine Mol_init_ions(this,nions,charges,positions)
      Class(T_Molecule),intent(inout) :: this
      integer                      :: nions
      real*8                       :: charges(:)
      real*8                       :: positions(:,:)

      if(this%has_ions()) &
         & call error("Mol_init_ions: ions already initialized")

      allocate(this%ions)
      this%ions=T_Ions(nions=nions)
      call this%ions%set_charges(charges)
      call this%ions%set_positions(positions)
   end subroutine

!Grid functions

   function Grid_has_weights(this) result(hasweights)
      Class(T_Grid),intent(in)     :: this
      Logical                      :: hasweights

      hasweights = allocated(this%weights)
   end function

   function Grid_has_cell(this) result(hascell)
      Class(T_Grid),intent(in)     :: this
      Logical                      :: hascell

      hascell = allocated(this%cell)
   end function

   subroutine Grid_set_weights(this,weights)
      Class(T_Grid),intent(inout) :: this
      real*8,intent(in)      :: weights(:)

      if(size(weights).ne.this%ngpt) &
         & call error("Grid_set_weights: size(weights) /= this%ngpt",size(weights),this%ngpt)

      if(.not.this%has_weights()) allocate(this%weights(this%ngpt))
      this%weights=weights(1:this%ngpt)
   end subroutine

   subroutine Grid_set_cell(this,cell)
      Class(T_Grid),intent(inout) :: this
      real*8,intent(in)      :: cell(:,:)

      if(size(cell).ne.9) &
         & call error("Grid_set_cell: shape(cell) incorrect",size(cell),9)

      if(.not.this%has_cell()) allocate(this%cell(3,3))
      this%cell=cell(1:3,1:3)
   end subroutine

!Ions functions
   function Ions_has_charges(this) result(hascharges)
      Class(T_Ions),intent(in)     :: this
      Logical                      :: hascharges

      hascharges = allocated(this%charges)
   end function

   function Ions_has_positions(this) result(haspositions)
      Class(T_Ions),intent(in)     :: this
      Logical                      :: haspositions

      haspositions = allocated(this%positions)
   end function

   subroutine Ions_set_charges(this,charges)
      Class(T_Ions),intent(inout) :: this
      real*8,intent(in)      :: charges(:)

      if(size(charges).ne.this%nions) &
         & call error("Ions_set_charges: size(charges) /= this%nions",size(charges),this%nions)

      if(.not.this%has_charges()) allocate(this%charges(this%nions))
      this%charges=charges(1:this%nions)
   end subroutine

   subroutine Ions_set_positions(this,positions)
      Class(T_Ions),intent(inout) :: this
      real*8,intent(in)      :: positions(:,:)

      if(size(positions).ne.3*this%nions) &
         & call error("Ions_set_positions: shape(positions incorrect",size(positions),3*this%nions)

      if(.not.this%has_positions()) allocate(this%positions(3,this%nions))
      this%positions=positions(1:3,1:this%nions)
   end subroutine

end
