!Types Module
!Contains all used types, including
!Molecule, Grid and Ions
module types_mod
   use output_mod
   implicit none
   private

   !Molecule Class, contains points, grid, density, nuclear potential and ions
   type public T_Molecule
      !#Grid Points
      integer                  :: ngpt
      !Grid, may contain weights and cell
      type(T_Grid),allocatable :: grid=T_Grid()
      !Electron Density
      real*8                   :: density(ngpt)
      !Nuclear Potential, alternative to ions
      real*8,allocatable       :: vnuc(:)
      !Ions,contains #ions, charges and positions
      type(T_Ions),allocatable :: ions=T_Ions()

    contains
      procedure :: set_grid      => Mol_set_grid
      procedure :: set_vnuc      => Mol_set_vnuc
      procedure :: set_ions      => Mol_set_ions

      procedure :: init_grid     => Mol_init_grid
      procedure :: init_ions     => Mol_init_ions

      procedure :: has_vnuc      => Mol_has_vnuc
      procedure :: has_ions      => Mol_has_ions
      procedure :: has_grid      => Mol_has_grid
   end type T_Molecule

   !Grid Class, contains points, weights and cell
   type public T_Grid
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
   type public T_Ions
      !#ions
      integer               :: nions
      !Charges
      real*8                :: charges(nions)
      !Positions,first index is x,y,z
      real*8                :: positions(3,nions)
   end type T_Ions


 contains

!Molecule functions

   function Mol_has_grid(this) result(hasgrid)
      Class(T_Molecule),intent(in) :: this
      Logical                      :: hasgrid
      
      hasgrid = allocated(this%grid) 
   end function

   function Mol_has_vnuc(this) result(hasvnuc)
      Class(T_Molecule),intent(in) :: this
      Logical                      :: hasvnuc
      
      hasvnuc = allocated(this%vnuc) 
   end function

   function Mol_has_ions(this) result(hasions)
      Class(T_Molecule),intent(in) :: this
      Logical                      :: hasions
      
      hasions = allocated(this%ions) 
   end function

   subroutine Mol_set_grid(this,grid)
      Class(T_Molecule),intent(in) :: this
      type(T_Grid),intent(in)      :: grid

      if(grid%ngpt.ne.this%ngpt) &
         & call error("Mol_set_grid "//str(this)//"grid%ngpt="//str(grid%ngpt)//"not equal ngpt="//str(ngpt))

      if(.not.this%has_grid) allocate(this%grid)
      this%grid=grid
   end subroutine

   subroutine Mol_set_vnuc(this,vnuc)
      Class(T_Molecule),intent(in) :: this
      real*8,intent(in)      :: vnuc(*)

      if(size(vnuc).ne.this%ngpt) &
         & call error("Mol_set_vnuc "//str(this)//"size(vnuc)="//str(size(vnuc))//"not equal ngpt="//str(ngpt))

      if(.not.this%has_vnuc) allocate(this%vnuc(this%ngpt))
      this%vnuc=vnuc
   end subroutine

   subroutine Mol_set_ions(this,ions)
      Class(T_Molecule),intent(in) :: this
      type(T_Ions),intent(in)      :: ions

      if(.not.this%has_ions) allocate(this%ions)
      this%ions=ions
   end subroutine

   subroutine Mol_init_grid(this,weights,cell)
      Class(T_Molecule),intent(in) :: this
      real*8,optional              :: weights(*)
      real*8,optional              :: cell(*)

      if(this%has_grid()) &
         & call error("Mol_init_grid "//str(this)//" grid already initialized")

      allocate(this%grid)
      this%grid=T_Grid(ngpt=this%ngpt)

      if(present(weights)) then
         this%grid%set_weights(weights)
      endif

      if(present(cell) then 
         this%grid%set_cell(cell)
      endif
   end subroutine

   subroutine Mol_init_ions(this,nions,charges,positions)
      Class(T_Molecule),intent(in) :: this
      integer                      :: nions
      real*8                       :: charges(*)
      real*8                       :: positions(*)

      if(this%has_ions()) &
         & call error("Mol_init_ions "//str(this)//" ions already initialized")

      allocate(this%ions)
      this%ions=T_Ions(nions=nions,charges=charges,positions=positions)
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
      Class(T_Grid),intent(in) :: this
      real*8,intent(in)      :: weights(*)

      if(size(weights).ne.this%ngpt) &
         & call error("Grid_set_weights "//str(this)//"size(weights)="//str(size(weights))//"not equal ngpt="//str(this%ngpt))

      if(.not.this%has_weights) allocate(this%weights(this%ngpt))
      this%weights=weights
   end subroutine

   subroutine Grid_set_cell(this,cell)
      Class(T_Grid),intent(in) :: this
      real*8,intent(in)      :: cell(*)

      if(shape(cell).ne.(/3,3/)) &
         & call error("Grid_set_cell "//str(this)//"shape(cell)="//str(shape(cell))//"incorrect")

      if(.not.this%has_cell) allocate(this%cell(3,3))
      this%cell=cell
   end subroutine

   

end
