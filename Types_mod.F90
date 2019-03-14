!! project: PEREMB
!!
!!
!! file: Types_mod.F90
!!
!! summary: Definition of all Data Types and their functions
!!
!! author: Christian Schwermann
!!
!! e-mail: c.schwermann@wwu.de
!!
!! date: 14/03/1019
!!
!! ALL RIGHTS RESERVED
!! Copyright © Christian Schwermann 2019
!!
module Types_mod
   use Precision_mod
   use Output_mod
   implicit none

   public :: grid_t, ions_t, molecule_t

   public :: Mol_has_grid, Mol_has_vnuc, Mol_has_ions, Mol_has_density, Mol_has_gradient
   public :: Mol_set_grid, Mol_set_vnuc, Mol_set_ions, Mol_set_spin, Mol_set_density, &
      & Mol_set_density_a, Mol_set_density_b, Mol_set_gradient, Mol_set_gradient_aa, &
      & Mol_set_gradient_ab, Mol_set_gradient_bb, Mol_init_grid, Mol_init_ions

   public :: Grid_has_weights, Grid_has_cell, Grid_has_positions
   public :: Grid_set_weights, Grid_set_cell, Grid_set_positions

   public :: Ions_has_charges, Ions_has_positions
   public :: Ions_set_charges, Ions_set_positions


   !!Grid type, contains points, weights and cell
   type :: grid_t
      !! #Grid Points
      integer               :: ngpt
      !! weights, optional
      real(kind=DP), allocatable    :: weights(:)
      !! cell vectors, optional
      real(kind=DP), allocatable    :: cell(:, :)
      !! @TODO need more cell information: zero point, #points in each vector direction @ENDTODO
      !! Positions,first index is x,y,z, optional
      real(kind=DP), allocatable    :: positions(:, :)
   end type grid_t

   !! Ions class, contains number, charges and positions
   type :: ions_t
      !! #ions
      integer               :: nions
      !! Charges
      real(kind=DP), allocatable    :: charges(:)
      !! Positions, first index is x,y,z
      real(kind=DP), allocatable    :: positions(:, :)
   end type ions_t

   !! Molecule type, contains points, grid, density, nuclear potential and ions
   type :: molecule_t
      !! #Grid Points
      integer                  :: ngpt
      !! Grid, may contain weights and cell
      type(grid_t), pointer     :: grid => null()
      !! Electron Density; Size=2*ngpt for spin polarized, beta is from ngpt+1:2*ngpt
      real(kind=DP), allocatable       :: density(:)
      !! Density Gradient; Size=3*ngpt for spin polarized, order is aa,ab,bb
      real(kind=DP), allocatable       :: gradient(:)
      !! Nuclear Potential, alternative to ions
      real(kind=DP), allocatable       :: vnuc(:)
      !! Ions, contains #ions, charges and positions
      type(ions_t), pointer     :: ions => null()
      !! Spin polarized or not
      logical                  :: spinpol = .false.
      !! Active or not
      logical                  :: active = .false.
   end type molecule_t

contains

!Molecule functions
   
   !! Return .true. if the molecule has a grid.
   pure function Mol_has_grid( this ) result( hasgrid )
      !! Input: molecule 
      type(molecule_t), intent(in) :: this
      !! Output: whether the molecule has a grid
      logical                      :: hasgrid

      hasgrid = Associated( this%grid )
   end function Mol_has_grid

   !! Return .true. if the molecule has a vnuc (nuclear potential).
   pure function Mol_has_vnuc( this ) result( hasvnuc )
      !! Input: molecule 
      type(molecule_t), intent(in) :: this
      !! Output: whether the molecule has a vnuc
      logical                      :: hasvnuc

      hasvnuc = Allocated( this%vnuc )
   end function Mol_has_vnuc

   !! Return .true. if the molecule has ions.
   pure function Mol_has_ions( this ) result( hasions )
      !! Input: molecule 
      type(molecule_t), intent(in) :: this
      !! Output: whether the molecule has ions
      logical                      :: hasions

      hasions = Associated( this%ions )
   end function Mol_has_ions

   !! Return .true. if the molecule has a density.
   pure function Mol_has_density( this ) result( hasdensity )
      !! Input: molecule 
      type(molecule_t), intent(in) :: this
      !! Output: whether the molecule has a density
      logical                      :: hasdensity

      hasdensity = Allocated( this%density )
   end function Mol_has_density

   !! Return .true. if the molecule has a density gradient.
   pure function Mol_has_gradient( this ) result( hasgradient )
      !! Input: molecule 
      type(molecule_t), intent(in) :: this
      !! Output: whether the molecule has a gradient
      logical                      :: hasgradient

      hasgradient = Allocated( this%gradient )
   end function Mol_has_gradient

   !! Set the grid of a molecule. Initialize a grid if none exists.
   subroutine Mol_set_grid( this, grid )
      !! Input: molecule
      type(molecule_t), intent(inout) :: this
      !! Input: grid, grid%ngpt has to be this%ngpt
      type(grid_t), intent(in)      :: grid

      if( grid%ngpt /= this%ngpt ) &
         & call Error( "Mol_set_grid: grid%ngpt /= this%ngpt", grid%ngpt, this%ngpt )

      if( .not. Mol_has_grid( this ) ) allocate( this%grid )
      this%grid = grid
   end subroutine Mol_set_grid

   !! Set the nuclear potential of a molecule. Initialize one if none exists.
   subroutine Mol_set_vnuc( this, vnuc )
      !! Input: molecule 
      type(molecule_t), intent(inout) :: this
      !! Input: nuclear potential, size has to be this%ngpt
      real(kind=DP), intent(in)      :: vnuc(:)
      !! Internal: start and end index of vnuc. 
      integer :: istart, iend

      if( Size( vnuc ) /= this%ngpt ) &
         & call Error( "Mol_set_vnuc: Size(vnuc) /= this%ngpt", Size( vnuc ), this%ngpt )

      istart = Lbound( vnuc, 1 )
      iend = Ubound( vnuc, 1 ) 

      if( .not. Mol_has_vnuc( this ) ) allocate( this%vnuc(1:this%ngpt) )
      this%vnuc(1:this%ngpt) = vnuc(istart:iend)
   end subroutine Mol_set_vnuc

   !! Set the ions of a molecule. Initialize if none exist.
   pure subroutine Mol_set_ions( this, ions )
      !! Input: molecule 
      type(molecule_t), intent(inout) :: this
      !! Input: ions 
      type(ions_t), intent(in)      :: ions

      if( .not. Mol_has_ions( this ) ) allocate( this%ions )
      this%ions = ions
   end subroutine Mol_set_ions

   !! Set the spin of a molecule. If density or gradient exist, they are deallocated!
   subroutine Mol_set_spin( this, spinpol )
      !! Input: molecule 
      type(molecule_t), intent(inout) :: this
      !! Input: spin ( .true. = spin polarized, .false. = unpolarized )
      logical, intent(in)             :: spinpol

      this%spinpol = spinpol

      if( Mol_has_density( this ) ) then 
         call Warning( "Mol_set_spin: setting spin after density! Deallocating density!" )
         deallocate( this%density )
      endif
      if( Mol_has_gradient( this ) ) then 
         call Warning( "Mol_set_spin: setting spin after gradient! Deallocating gradient!" )
         deallocate( this%gradient )
      endif
   end subroutine Mol_set_spin

   subroutine Mol_set_density( this, density )
      type(molecule_t), intent(inout) :: this
      real(kind=DP), intent(in)            :: density(:)
      integer                      :: dsize, istart, iend

      if( this%spinpol ) then
         dsize = 2 * this%ngpt
      else
         dsize = this%ngpt
      end if

      if( Size( density ) /= dsize ) &
         & call Error( "Mol_set_density: Size(density) /= 2*this%ngpt", Size(density), dsize )

      istart = Lbound( density, 1 )
      iend = Ubound( density, 1 )

      if( .not. Mol_has_density( this ) ) allocate( this%density(1:dsize) )
      this%density(1:dsize) = density(istart:iend)
   end subroutine Mol_set_density

   subroutine Mol_set_density_a( this, density )
      type(molecule_t), intent(inout) :: this
      real(kind=DP), intent(in)            :: density(:)
      integer                        :: dsize, istart, iend

      dsize = 2 * this%ngpt

      if( .not. this%spinpol ) &
         & call Error( "Mol_set_density_a: trying to set alpha density for non spin polarized molecule!" )

      if( Size( density ) /= this%ngpt ) &
         & call Error( "Mol_set_density_a: Size(density) /= this%ngpt", Size( density ), this%ngpt )

      istart = Lbound( density, 1 )
      iend = Ubound( density, 1 )

      if( .not. Mol_has_density( this ) ) allocate( this%density(1:dsize) )
      this%density(1:this%ngpt) = density(istart:iend)
   end subroutine Mol_set_density_a

   subroutine Mol_set_density_b( this, density )
      type(molecule_t), intent(inout) :: this
      real(kind=DP), intent(in)            :: density(:)
      integer                        :: dsize, istart, iend

      dsize = 2 * this%ngpt

      if( .not. this%spinpol ) &
         & call Error( "Mol_set_density_b: trying to set beta density for non spin polarized molecule!" )

      if( Size( density ) /= this%ngpt ) &
         & call Error( "Mol_set_density_b: Size(density) /= this%ngpt", Size( density ), this%ngpt )

      istart = Lbound( density, 1 )
      iend = Ubound( density, 1 )

      if( .not. Mol_has_density( this ) ) allocate( this%density(1:dsize) )
      this%density(this%ngpt+1:dsize) = density(istart:iend)
   end subroutine Mol_set_density_b

   subroutine Mol_set_gradient( this, gradient )
      type(molecule_t), intent(inout) :: this
      real(kind=DP), intent(in)            :: gradient(:)
      integer                      :: dsize, istart, iend

      if( this%spinpol ) then
         dsize = 3 * this%ngpt
      else
         dsize = this%ngpt
      end if

      if( Size( gradient ) /= dsize ) &
         & call Error( "Mol_set_gradient: Size(gradient) /= 3*this%ngpt", Size( gradient ), dsize )

      istart = Lbound( gradient, 1 )
      iend = Ubound( gradient, 1 )

      if( .not. Mol_has_gradient( this ) ) allocate( this%gradient(dsize) )
      this%gradient(1:dsize) = gradient(istart:iend)
   end subroutine Mol_set_gradient

   subroutine Mol_set_gradient_aa( this, gradient )
      type(molecule_t), intent(inout) :: this
      real(kind=DP), intent(in)            :: gradient(:)
      integer                      :: dsize, istart, iend

      dsize = 3 * this%ngpt

      if( .not. this%spinpol ) &
         & call Error( "Mol_set_gradient_aa: trying to set alpha gradient for non spin polarized molecule!" )

      if( Size( gradient ) /= dsize ) &
         & call Error( "Mol_set_gradient_aa: Size(gradient) /= this%ngpt", Size( gradient ), this%ngpt )

      istart = Lbound( gradient, 1 )
      iend = Ubound( gradient, 1 )

      if( .not. Mol_has_gradient( this ) ) allocate( this%gradient(dsize) )
      this%gradient(1:this%ngpt) = gradient(istart:iend)
   end subroutine Mol_set_gradient_aa

   subroutine Mol_set_gradient_ab( this, gradient )
      type(molecule_t), intent(inout) :: this
      real(kind=DP), intent(in)            :: gradient(:)
      integer                      :: dsize, istart, iend

      dsize= 3 * this%ngpt

      if( .not. this%spinpol ) &
         & call Error( "Mol_set_gradient_ab: trying to set alpha*beta gradient for non spin polarized molecule!" )

      if( Size( gradient ) /= this%ngpt ) &
         & call Error( "Mol_set_gradient_ab: Size(gradient) /= this%ngpt", Size( gradient ), this%ngpt )

      istart = Lbound( gradient, 1 )
      iend = Ubound( gradient, 1 )

      if( .not. Mol_has_gradient( this ) ) allocate( this%gradient(dsize) )
      this%gradient(this%ngpt+1:2*this%ngpt) = gradient(istart:iend)
   end subroutine Mol_set_gradient_ab

   subroutine Mol_set_gradient_bb( this, gradient )
      type(molecule_t), intent(inout) :: this
      real(kind=DP), intent(in)            :: gradient(:)
      integer                      :: dsize, istart, iend

      dsize = 3 * this%ngpt

      if( .not. this%spinpol ) &
         & call Error( "Mol_set_gradient_bb: trying to set beta gradient for non spin polarized molecule!" )

      if( Size( gradient ) /= this%ngpt ) &
         & call Error( "Mol_set_gradient_bb: Size(gradient) /= this%ngpt", Size( gradient ), this%ngpt )

      istart = Lbound( gradient, 1 )
      iend = Ubound( gradient, 1 )

      if( .not. Mol_has_gradient( this ) ) Allocate( this%gradient(dsize) )
      this%gradient(2*this%ngpt+1:dsize) = gradient(istart:iend)
   end subroutine Mol_set_gradient_bb

   subroutine Mol_init_grid( this, weights, cell, positions )
      type(molecule_t), intent(inout) :: this
      real(kind=DP), optional, intent(in)  :: weights(:)
      real(kind=DP), optional, intent(in)  :: cell(:, :)
      real(kind=DP), optional, intent(in)  :: positions(:, :)

      if( Mol_has_grid( this ) ) &
         & call Error( "Mol_init_grid: grid already initialized" )

      allocate( this%grid )
      this%grid = grid_t( ngpt = this%ngpt )

      if( Present( weights ) ) then
         call Grid_set_weights( this%grid, weights )
      end if

      if( Present( cell ) ) then
         call Grid_set_cell( this%grid, cell )
      end if

      if( Present( positions ) ) then
         call Grid_set_positions( this%grid, positions )
      end if
   end subroutine Mol_init_grid

   subroutine Mol_init_ions( this, nions, charges, positions )
      type(molecule_t), intent(inout) :: this
      integer, intent(in)          :: nions
      real(kind=DP), intent(in)           :: charges(:)
      real(kind=DP), intent(in)           :: positions(:, :)

      if( Mol_has_ions( this ) ) &
         & call Error( "Mol_init_ions: ions already initialized" )

      allocate( this%ions )
      this%ions = ions_t( nions = nions )
      call Ions_set_charges( this%ions, charges )
      call Ions_set_positions( this%ions, positions )
   end subroutine Mol_init_ions

!Grid functions

   pure function Grid_has_weights( this ) result( hasweights )
      type(grid_t), intent(in)     :: this
      logical                      :: hasweights

      hasweights = Allocated( this%weights )
   end function Grid_has_weights

   pure function Grid_has_cell( this ) result( hascell )
      type(grid_t), intent(in)     :: this
      logical                      :: hascell

      hascell = Allocated( this%cell )
   end function Grid_has_cell

   pure function Grid_has_positions( this ) result( haspositions )
      type(grid_t), intent(in)     :: this
      logical                      :: haspositions

      haspositions = Allocated( this%positions )
   end function Grid_has_positions

   subroutine Grid_set_weights( this, weights )
      type(grid_t), intent(inout) :: this
      real(kind=DP), intent(in)      :: weights(:)
      integer :: istart, iend

      if( Size( weights ) /= this%ngpt ) &
         & call Error( "Grid_set_weights: Size(weights) /= this%ngpt", Size( weights ), this%ngpt )

      istart = Lbound( weights, 1 )
      iend = Ubound( weights, 1 )

      if( .not. Grid_has_weights( this ) ) allocate( this%weights(1:this%ngpt) )
      this%weights(1:this%ngpt) = weights(istart:iend)
   end subroutine Grid_set_weights

   subroutine Grid_set_cell( this, cell )
      type(grid_t), intent(inout) :: this
      real(kind=DP), intent(in)      :: cell(:, :)
      integer :: istart1,iend1,istart2,iend2

      if( Size( cell ) /= 9 ) &
         & call Error( "Grid_set_cell: shape(cell) incorrect", Size( cell ), 9 )

      istart1 = Lbound( cell, 1 )
      iend1 = Ubound( cell, 1 )
      istart2 = Lbound( cell, 2 )
      iend2 = Ubound( cell, 2 )

      if( .not. Grid_has_cell( this ) ) allocate( this%cell(1:3, 1:3) )
      this%cell(1:3, 1:3) = cell(istart1:iend1, istart2:iend2)
   end subroutine Grid_set_cell

   subroutine Grid_set_positions( this, positions )
      type(grid_t), intent(inout) :: this
      real(kind=DP), intent(in)      :: positions(:, :)
      integer :: istart1, iend1, istart2, iend2

      if( Size( positions ) /= 3 * this%ngpt ) &
         & call Error( "Grid_set_positions: shape(positions) incorrect", Size( positions ), 3 * this%ngpt )

      istart1 = Lbound( positions, 1 )
      iend1 = Ubound( positions, 1 )
      istart2 = Lbound( positions, 2 )
      iend2 = Ubound( positions, 2 )

      if( .not. Grid_has_positions( this ) ) allocate( this%positions(1:3, 1:this%ngpt) )
      this%positions(1:3, 1:this%ngpt) = positions(istart1:iend1, istart2:iend2)
   end subroutine Grid_set_positions

!Ions functions
   pure function Ions_has_charges( this ) result( hascharges )
      type(ions_t), intent(in)     :: this
      logical                      :: hascharges

      hascharges = Allocated( this%charges )
   end function Ions_has_charges

   pure function Ions_has_positions( this ) result( haspositions )
      type(ions_t), intent(in)     :: this
      logical                      :: haspositions

      haspositions = Allocated( this%positions )
   end function Ions_has_positions

   subroutine Ions_set_charges( this, charges )
      type(ions_t), intent(inout) :: this
      real(kind=DP), intent(in)      :: charges(:)
      integer :: istart, iend

      if( Size( charges ) /= this%nions ) &
         & call Error("Ions_set_charges: Size(charges) /= this%nions", Size(charges), this%nions )

      istart = Lbound( charges, 1 )
      iend = Ubound( charges, 1 )

      if( .not. Ions_has_charges( this ) ) allocate( this%charges(1:this%nions) )
      this%charges(1:this%nions) = charges(istart:iend)
   end subroutine Ions_set_charges

   subroutine Ions_set_positions( this, positions )
      type(ions_t), intent(inout) :: this
      real(kind=DP), intent(in)      :: positions(:, :)
      integer :: istart1, iend1, istart2, iend2

      if( Size( positions ) /= 3 * this%nions ) &
         & call Error( "Ions_set_positions: shape(positions incorrect", Size( positions ), 3 * this%nions )

      istart1 = Lbound( positions, 1 )
      iend1 = Ubound( positions, 1 )
      istart2 = Lbound( positions, 2 )
      iend2 = Ubound( positions, 2 )

      if( .not. Ions_has_positions( this ) ) allocate( this%positions(1:3, 1:this%nions) )
      this%positions(1:3, 1:this%nions) = positions(istart1:iend1, istart2:iend2)
   end subroutine Ions_set_positions

end module Types_mod
