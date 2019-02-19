!Types Module
!Contains all used types, including
!Molecule, Grid and Ions
module types_mod
   use output_mod
   implicit none
   public

   !Grid type, contains points, weights and cell
   type, public :: T_Grid
      !#Grid Points
      integer               :: ngpt
      !weights,optional
      real*8,allocatable    :: weights(:)
      !cell vectors,optional
      real*8,allocatable    :: cell(:,:)
   end type T_Grid

   !Ions class, contains number, charges and positions
   type, public :: T_Ions
      !#ions
      integer               :: nions
      !Charges
      real*8,allocatable    :: charges(:)
      !Positions,first index is x,y,z
      real*8,allocatable    :: positions(:,:)
   end type T_Ions

   !Molecule type, contains points, grid, density, nuclear potential and ions
   type, public :: T_Molecule
      !#Grid Points
      integer                  :: ngpt
      !Grid, may contain weights and cell
      type(T_Grid),pointer     :: grid => null()
      !Electron Density; size=2*ngpt for spin polarized, beta is from ngpt+1:2*ngpt
      real*8,allocatable       :: density(:)
      !Nuclear Potential, alternative to ions
      real*8,allocatable       :: vnuc(:)
      !Ions,contains #ions, charges and positions
      type(T_Ions),pointer     :: ions => null()
      !Spin polarized or not
      logical                  :: spinpol=.false.
   end type T_Molecule

 contains

!Molecule functions

   function Mol_has_grid(this) result(hasgrid)
      type(T_Molecule),intent(in) :: this
      Logical                      :: hasgrid

      hasgrid = associated(this%grid)
   end function

   function Mol_has_vnuc(this) result(hasvnuc)
      type(T_Molecule),intent(in) :: this
      Logical                      :: hasvnuc

      hasvnuc = allocated(this%vnuc)
   end function

   function Mol_has_ions(this) result(hasions)
      type(T_Molecule),intent(in) :: this
      Logical                      :: hasions

      hasions = associated(this%ions)
   end function

   function Mol_has_density(this) result(hasdensity)
      type(T_Molecule),intent(in) :: this
      Logical                      :: hasdensity

      hasdensity = allocated(this%density)
   end function

   subroutine Mol_set_grid(this,grid)
      type(T_Molecule),intent(inout) :: this
      type(T_Grid),intent(in)      :: grid

      if(grid%ngpt.ne.this%ngpt) &
         & call error("Mol_set_grid: grid%ngpt /= this%ngpt",grid%ngpt,this%ngpt)

      if(.not.Mol_has_grid(this)) allocate(this%grid)
      this%grid=grid
   end subroutine

   subroutine Mol_set_vnuc(this,vnuc)
      type(T_Molecule),intent(inout) :: this
      real*8,intent(in)      :: vnuc(:)

      if(size(vnuc).ne.this%ngpt) &
         & call error("Mol_set_vnuc: size(vnuc) /= this%ngpt",size(vnuc),this%ngpt)

      if(.not.Mol_has_vnuc(this)) allocate(this%vnuc(this%ngpt))
      this%vnuc=vnuc(1:this%ngpt)
   end subroutine

   subroutine Mol_set_ions(this,ions)
      type(T_Molecule),intent(inout) :: this
      type(T_Ions),intent(in)      :: ions

      if(.not.Mol_has_ions(this)) allocate(this%ions)
      this%ions=ions
   end subroutine

   subroutine Mol_set_density(this,density)
      type(T_Molecule),intent(inout) :: this
      real*8,intent(in)            :: density(:)
      integer                      :: dsize

      if(this%spinpol) then
         dsize=2*this%ngpt
      else
         dsize=this%ngpt
      endif

      if(size(density).ne.dsize) &
         & call error("Mol_set_density: size(density) /= spins*this%ngpt",size(density),dsize)

      if(.not.Mol_has_density(this)) allocate(this%density(dsize))
      this%density=density(1:dsize)
   end subroutine

   subroutine Mol_set_density_a(this,density)
      type(T_Molecule),intent(inout) :: this
      real*8,intent(in)            :: density(:)
      integer                        :: dsize

      dsize=2*this%ngpt

      if(.not.this%spinpol) &
         & call error("Mol_set_density_a: trying to set alpha density for non spin polarized molecule!")

      if(size(density).ne.this%ngpt) &
         & call error("Mol_set_density_a: size(density) /= this%ngpt",size(density),this%ngpt)

      if(.not.Mol_has_density(this)) allocate(this%density(dsize))
      this%density(1:this%ngpt)=density(1:this%ngpt)
   end subroutine

   subroutine Mol_set_density_b(this,density)
      type(T_Molecule),intent(inout) :: this
      real*8,intent(in)            :: density(:)
      integer                        :: dsize

      dsize=2*this%ngpt

      if(.not.this%spinpol) &
         & call error("Mol_set_density_b: trying to set beta density for non spin polarized molecule!")

      if(size(density).ne.this%ngpt) &
         & call error("Mol_set_density_b: size(density) /= this%ngpt",size(density),this%ngpt)

      if(.not.Mol_has_density(this)) allocate(this%density(dsize))
      this%density(this%ngpt+1:dsize)=density(1:this%ngpt)
   end subroutine

   subroutine Mol_init_grid(this,weights,cell)
      type(T_Molecule),intent(inout) :: this
      real*8,optional              :: weights(:)
      real*8,optional              :: cell(:,:)

      if(Mol_has_grid(this)) &
         & call error("Mol_init_grid: grid already initialized")

      allocate(this%grid)
      this%grid=T_Grid(ngpt=this%ngpt)

      if(present(weights)) then
         call Grid_set_weights(this%grid,weights)
      endif

      if(present(cell)) then
         call Grid_set_cell(this%grid,cell)
      endif
   end subroutine

   subroutine Mol_init_ions(this,nions,charges,positions)
      type(T_Molecule),intent(inout) :: this
      integer                      :: nions
      real*8                       :: charges(:)
      real*8                       :: positions(:,:)

      if(Mol_has_ions(this)) &
         & call error("Mol_init_ions: ions already initialized")

      allocate(this%ions)
      this%ions=T_Ions(nions=nions)
      call Ions_set_charges(this%ions,charges)
      call Ions_set_positions(this%ions,positions)
   end subroutine

!Grid functions

   function Grid_has_weights(this) result(hasweights)
      type(T_Grid),intent(in)     :: this
      Logical                      :: hasweights

      hasweights = allocated(this%weights)
   end function

   function Grid_has_cell(this) result(hascell)
      type(T_Grid),intent(in)     :: this
      Logical                      :: hascell

      hascell = allocated(this%cell)
   end function

   subroutine Grid_set_weights(this,weights)
      type(T_Grid),intent(inout) :: this
      real*8,intent(in)      :: weights(:)

      if(size(weights).ne.this%ngpt) &
         & call error("Grid_set_weights: size(weights) /= this%ngpt",size(weights),this%ngpt)

      if(.not.Grid_has_weights(this)) allocate(this%weights(this%ngpt))
      this%weights=weights(1:this%ngpt)
   end subroutine

   subroutine Grid_set_cell(this,cell)
      type(T_Grid),intent(inout) :: this
      real*8,intent(in)      :: cell(:,:)

      if(size(cell).ne.9) &
         & call error("Grid_set_cell: shape(cell) incorrect",size(cell),9)

      if(.not.Grid_has_cell(this)) allocate(this%cell(3,3))
      this%cell=cell(1:3,1:3)
   end subroutine

!Ions functions
   function Ions_has_charges(this) result(hascharges)
      type(T_Ions),intent(in)     :: this
      Logical                      :: hascharges

      hascharges = allocated(this%charges)
   end function

   function Ions_has_positions(this) result(haspositions)
      type(T_Ions),intent(in)     :: this
      Logical                      :: haspositions

      haspositions = allocated(this%positions)
   end function

   subroutine Ions_set_charges(this,charges)
      type(T_Ions),intent(inout) :: this
      real*8,intent(in)      :: charges(:)

      if(size(charges).ne.this%nions) &
         & call error("Ions_set_charges: size(charges) /= this%nions",size(charges),this%nions)

      if(.not.Ions_has_charges(this)) allocate(this%charges(this%nions))
      this%charges=charges(1:this%nions)
   end subroutine

   subroutine Ions_set_positions(this,positions)
      type(T_Ions),intent(inout) :: this
      real*8,intent(in)      :: positions(:,:)

      if(size(positions).ne.3*this%nions) &
         & call error("Ions_set_positions: shape(positions incorrect",size(positions),3*this%nions)

      if(.not.Ions_has_positions(this)) allocate(this%positions(3,this%nions))
      this%positions=positions(1:3,1:this%nions)
   end subroutine

end
