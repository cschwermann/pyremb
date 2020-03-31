!! Author:  Christian Schwermann
!! E-mail:  c.schwermann@wwu.de
!! Date:    14/03/2019
!! Project: PEREMB 
!! File:    Types_mod.F90 
!! Copyright: © 2019 Christian Schwermann, ALL RIGHTS RESERVED
!!
!!***********************************************************************
!!
!! Definition of all Data Types and their functions
!!
!!***********************************************************************
module Types_mod
   use Precision_mod ! IEEE precision
   use Output_mod    ! Error and warning
   implicit none

   ! Grid type, Ions type, Molecule type
   public :: grid_t, ions_t, molecule_t

   ! Molecule functions and subroutines
   public :: Mol_has_grid, Mol_has_vnuc, Mol_has_ions, Mol_has_density, Mol_has_gradient
   public :: Mol_set_grid, Mol_set_vnuc, Mol_set_ions, Mol_set_spin, Mol_set_density, &
      & Mol_set_density_a, Mol_set_density_b, Mol_set_gradient, Mol_set_gradient_aa, &
      & Mol_set_gradient_ab, Mol_set_gradient_bb, Mol_init_grid, Mol_init_ions

   ! Grid functions and subroutines
   public :: Grid_has_weights, Grid_has_cell, Grid_has_positions
   public :: Grid_set_weights, Grid_set_cell, Grid_set_positions

   ! Ions functions and subroutines
   public :: Ions_has_charges, Ions_has_positions
   public :: Ions_set_charges, Ions_set_positions

   !!********************************************************************
   !! Grid type, contains points, weights, cell and positions
   !!********************************************************************
   type :: grid_t
      !! Number of grid points
      integer                       :: ngpt
      !! Weights; size = ngpt, optional
      real(kind=DP), allocatable    :: weights(:)
      !! Cell vectors; size = (3, 3), optional
      real(kind=DP), allocatable    :: cell(:, :)
      !TODO need more cell information: zero point, #points in each vector direction
      !! Positions; size = (3, ngpt), first index is x,y,z, optional
      real(kind=DP), allocatable    :: positions(:, :)
      !!

   end type grid_t


   !!********************************************************************
   !! Ions type, contains number, charges and positions
   !!********************************************************************
   type :: ions_t
      !! Number of ions
      integer                       :: nions
      !! Charges; size = nions
      real(kind=DP), allocatable    :: charges(:)
      !! Positions; size = (3, nions), first index is x,y,z
      real(kind=DP), allocatable    :: positions(:, :)
      !!

   end type ions_t


   !!********************************************************************
   !! Molecule type, contains points, grid, density, nuclear potential and ions
   !!********************************************************************
   type :: molecule_t
      !! Number of grid points
      integer                       :: ngpt
      !! Grid, may contain weights, cell and positions
      type(grid_t), pointer         :: grid => null()
      !! Electron density; size = ngpt for unpolarized or 2 * ngpt for spin polarized, beta is from ngpt + 1:2 * ngpt
      real(kind=DP), allocatable    :: density(:)
      !! Density gradient \(|\nabla\rho|^2\); size = npgt for unpolarized or 3 * ngpt for spin polarized, order is aa,ab,bb
      real(kind=DP), allocatable    :: gradient(:)
      !! Nuclear potential; size = ngpt, alternative to ions
      real(kind=DP), allocatable    :: vnuc(:)
      !! Ions, contains number of ions, charges and positions
      type(ions_t), pointer         :: ions => null()
      !! Spin polarized or not
      logical                       :: spinpol = .false.
      !! Active or not
      logical                       :: active = .false.
      !!

   end type molecule_t


contains


   !********************************************************************
   ! Molecule functions
   !********************************************************************
   
   !! Return .true. if the molecule has a grid.
   pure function Mol_has_grid( this ) result( hasgrid )
      !! Input: molecule 
      type(molecule_t), intent(in) :: this
      !! Output: whether the molecule has a grid
      logical                      :: hasgrid
      !!

      hasgrid = Associated( this%grid )

   end function Mol_has_grid


   !! Return .true. if the molecule has a vnuc (nuclear potential).
   pure function Mol_has_vnuc( this ) result( hasvnuc )
      !! Input: molecule 
      type(molecule_t), intent(in) :: this
      !! Output: whether the molecule has a vnuc
      logical                      :: hasvnuc
      !!

      hasvnuc = Allocated( this%vnuc )

   end function Mol_has_vnuc


   !! Return .true. if the molecule has ions.
   pure function Mol_has_ions( this ) result( hasions )
      !! Input: molecule 
      type(molecule_t), intent(in) :: this
      !! Output: whether the molecule has ions
      logical                      :: hasions
      !!

      hasions = Associated( this%ions )

   end function Mol_has_ions


   !! Return .true. if the molecule has a density.
   pure function Mol_has_density( this ) result( hasdensity )
      !! Input: molecule 
      type(molecule_t), intent(in) :: this
      !! Output: whether the molecule has a density
      logical                      :: hasdensity
      !!

      hasdensity = Allocated( this%density )

   end function Mol_has_density


   !! Return .true. if the molecule has a density gradient.
   pure function Mol_has_gradient( this ) result( hasgradient )
      !! Input: molecule 
      type(molecule_t), intent(in) :: this
      !! Output: whether the molecule has a gradient
      logical                      :: hasgradient
      !!

      hasgradient = Allocated( this%gradient )

   end function Mol_has_gradient


   !! Set the grid of a molecule. Initialize a grid if none exists.
   subroutine Mol_set_grid( this, grid )
      !! Input: molecule
      type(molecule_t), intent(inout) :: this
      !! Input: grid, ngpt has to be ngpt of molecule
      type(grid_t), intent(in)        :: grid
      !!

      if( grid%ngpt /= this%ngpt ) &
         & call Error( "Mol_set_grid: grid%ngpt /= this%ngpt", grid%ngpt, this%ngpt )

      if( .not. Mol_has_grid( this ) ) allocate( this%grid )
      this%grid = grid

   end subroutine Mol_set_grid


   !! Set the nuclear potential of a molecule. Initialize one if none exists.
   subroutine Mol_set_vnuc( this, vnuc )
      !! Input: molecule 
      type(molecule_t), intent(inout) :: this
      !! Input: nuclear potential, 1D array of size ngpt
      real(kind=DP), intent(in)       :: vnuc(:)
      !! Internal: start and end index of vnuc
      integer :: istart, iend
      !!

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
      type(ions_t), intent(in)        :: ions
      !!

      if( .not. Mol_has_ions( this ) ) allocate( this%ions )
      this%ions = ions

   end subroutine Mol_set_ions


   !! Set the spin of a molecule. If density or gradient exist, they are deallocated!
   subroutine Mol_set_spin( this, spinpol )
      !! Input: molecule 
      type(molecule_t), intent(inout) :: this
      !! Input: spin ( .true. = spin polarized, .false. = unpolarized )
      logical, intent(in)             :: spinpol
      !!

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


   !! Set the density of a molecule. Initialize density if none exists.
   subroutine Mol_set_density( this, density )
      !! Input: molecule 
      type(molecule_t), intent(inout) :: this
      !! Input: density, 1D array of length ngpt (or 2*ngpt if spin polarized; order is a, b) 
      real(kind=DP), intent(in)       :: density(:)
      !! Internal :: supposed size of density, start and end of input density
      integer                         :: dsize, istart, iend
      !!

      if( this%spinpol ) then
         dsize = 2 * this%ngpt
      else
         dsize = this%ngpt
      end if

      if( Size( density ) /= dsize ) &
         & call Error( "Mol_set_density: Size(density) /= 2 * this%ngpt", Size(density), dsize )

      istart = Lbound( density, 1 )
      iend = Ubound( density, 1 )

      if( .not. Mol_has_density( this ) ) allocate( this%density(1:dsize) )
      this%density(1:dsize) = density(istart:iend)

   end subroutine Mol_set_density


   !! Set the alpha spin density of a molecule. Initialize density if none exists.
   subroutine Mol_set_density_a( this, density )
      !! Input: molecule 
      type(molecule_t), intent(inout) :: this
      !! Input: density, 1D array of length ngpt 
      real(kind=DP), intent(in)       :: density(:)
      !! Internal :: supposed size of density, start and end of input density
      integer                         :: dsize, istart, iend
      !!

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


   !! Set the beta spin density of a molecule. Initialize density if none exists.
   subroutine Mol_set_density_b( this, density )
      !! Input: molecule 
      type(molecule_t), intent(inout) :: this
      !! Input: density, 1D array of length ngpt 
      real(kind=DP), intent(in)       :: density(:)
      !! Internal :: supposed size of density, start and end of input density
      integer                         :: dsize, istart, iend
      !!

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


   !! Set the density gradient of a molecule. Initialize gradient if none exists.
   subroutine Mol_set_gradient( this, gradient )
      !! Input: molecule 
      type(molecule_t), intent(inout) :: this
      !! Input: gradient, 1D array of length ngpt (3* ngpt if spin polarized; order is aa, ab, bb)
      real(kind=DP), intent(in)       :: gradient(:)
      !! Internal :: supposed size of gradient, start and end of input gradient
      integer                         :: dsize, istart, iend
      !!

      if( this%spinpol ) then
         dsize = 3 * this%ngpt
      else
         dsize = this%ngpt
      end if

      if( Size( gradient ) /= dsize ) &
         & call Error( "Mol_set_gradient: Size(gradient) /= 3 * this%ngpt", Size( gradient ), dsize )

      istart = Lbound( gradient, 1 )
      iend = Ubound( gradient, 1 )

      if( .not. Mol_has_gradient( this ) ) allocate( this%gradient(dsize) )
      this%gradient(1:dsize) = gradient(istart:iend)

   end subroutine Mol_set_gradient


   !! Set the aa density gradient of a molecule. Initialize gradient if none exists.
   subroutine Mol_set_gradient_aa( this, gradient )
      !! Input: molecule 
      type(molecule_t), intent(inout) :: this
      !! Input: gradient, 1D array of length ngpt
      real(kind=DP), intent(in)       :: gradient(:)
      !! Internal :: supposed size of gradient, start and end of input gradient
      integer                         :: dsize, istart, iend
      !!

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


   !! Set the aa density gradient of a molecule. Initialize gradient if none exists.
   subroutine Mol_set_gradient_ab( this, gradient )
      !! Input: molecule 
      type(molecule_t), intent(inout) :: this
      !! Input: gradient, 1D array of length ngpt
      real(kind=DP), intent(in)       :: gradient(:)
      !! Internal :: supposed size of gradient, start and end of input gradient
      integer                         :: dsize, istart, iend
      !!

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


   !! Set the aa density gradient of a molecule. Initialize gradient if none exists.
   subroutine Mol_set_gradient_bb( this, gradient )
      !! Input: molecule 
      type(molecule_t), intent(inout) :: this
      !! Input: gradient, 1D array of length ngpt
      real(kind=DP), intent(in)       :: gradient(:)
      !! Internal :: supposed size of gradient, start and end of input gradient
      integer                         :: dsize, istart, iend
      !!

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


   !! Initialize the grid of a molecule. This requires nothing, but accepts weights, cell information and grid positions.
   subroutine Mol_init_grid( this, weights, cell, positions )
      !! Input: molecule 
      type(molecule_t), intent(inout)     :: this
      !! Input: grid weights, size ngpt
      real(kind=DP), optional, intent(in) :: weights(:)
      !! Input: unit cell, size (3, 3)
      real(kind=DP), optional, intent(in) :: cell(:, :)
      !! Input: grid weights, size (3, ngpt); first index is x, y, z
      real(kind=DP), optional, intent(in) :: positions(:, :)
      !!

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


   !! Initialize the ions of a molecule. This requires #ions, their charges and xyz positions.
   subroutine Mol_init_ions( this, nions, charges, positions )
      !! Input: molecule 
      type(molecule_t), intent(inout) :: this
      !! Input: number of ions
      integer, intent(in)             :: nions
      !! Input: nuclear (effective) charges, size nions
      real(kind=DP), intent(in)       :: charges(:)
      !! Input: nuclear positions, size (3, nions); first index is x, y, z
      real(kind=DP), intent(in)       :: positions(:, :)
      !!

      if( Mol_has_ions( this ) ) &
         & call Error( "Mol_init_ions: ions already initialized" )

      allocate( this%ions )
      this%ions = ions_t( nions = nions )
      call Ions_set_charges( this%ions, charges )
      call Ions_set_positions( this%ions, positions )

   end subroutine Mol_init_ions


   !********************************************************************
   ! Grid functions
   !********************************************************************

   !! Return .true. if the grid has weights.
   pure function Grid_has_weights( this ) result( hasweights )
      !! Input: grid
      type(grid_t), intent(in) :: this
      !! Output: whether the grid has weights
      logical                  :: hasweights
      !!

      hasweights = Allocated( this%weights )

   end function Grid_has_weights


   !! Return .true. if the grid has a cell.
   pure function Grid_has_cell( this ) result( hascell )
      !! Input: grid
      type(grid_t), intent(in) :: this
      !! Output: whether the grid has a cell
      logical                  :: hascell
      !!

      hascell = Allocated( this%cell )

   end function Grid_has_cell


   !! Return .true. if the grid has positions.
   pure function Grid_has_positions( this ) result( haspositions )
      !! Input: grid
      type(grid_t), intent(in) :: this
      !! Output: whether the grid has a cell
      logical                  :: haspositions
      !!

      haspositions = Allocated( this%positions )

   end function Grid_has_positions


   !! Set the weights of a grid. Initialize weights if none exist.
   subroutine Grid_set_weights( this, weights )
      !! Input: grid
      type(grid_t), intent(inout) :: this
      !! Input: weights, 1D array of size ngpt
      real(kind=DP), intent(in)   :: weights(:)
      !! Internal :: start and end of input weights
      integer                     :: istart, iend
      !!

      if( Size( weights ) /= this%ngpt ) &
         & call Error( "Grid_set_weights: Size(weights) /= this%ngpt", Size( weights ), this%ngpt )

      istart = Lbound( weights, 1 )
      iend = Ubound( weights, 1 )

      if( .not. Grid_has_weights( this ) ) allocate( this%weights(1:this%ngpt) )
      this%weights(1:this%ngpt) = weights(istart:iend)

   end subroutine Grid_set_weights


   !! Set the cell of a grid. Initialize cell if none exists.
   subroutine Grid_set_cell( this, cell )
      !! Input: grid
      type(grid_t), intent(inout) :: this
      !! Input: cell, 2D array of size (3, 3), contains cell vectors
      real(kind=DP), intent(in)   :: cell(:, :)
      !! Internal :: start and end of input cell for each dimension
      integer                     :: istart1, iend1, istart2, iend2
      !!

      if( Size( cell ) /= 9 ) &
         & call Error( "Grid_set_cell: shape(cell) incorrect", Size( cell ), 9 )

      istart1 = Lbound( cell, 1 )
      iend1 = Ubound( cell, 1 )
      istart2 = Lbound( cell, 2 )
      iend2 = Ubound( cell, 2 )

      if( .not. Grid_has_cell( this ) ) allocate( this%cell(1:3, 1:3) )
      this%cell(1:3, 1:3) = cell(istart1:iend1, istart2:iend2)

   end subroutine Grid_set_cell


   !! Set the positions of a grid. Initialize positions if none exist.
   subroutine Grid_set_positions( this, positions )
      !! Input: grid
      type(grid_t), intent(inout) :: this
      !! Input: cell, 2D array of size (3, ngpt), first index is x, y, z
      real(kind=DP), intent(in)   :: positions(:, :)
      !! Internal :: start and end of input positions for each dimension
      integer                     :: istart1, iend1, istart2, iend2
      !!

      if( Size( positions ) /= 3 * this%ngpt ) &
         & call Error( "Grid_set_positions: shape(positions) incorrect", Size( positions ), 3 * this%ngpt )

      istart1 = Lbound( positions, 1 )
      iend1 = Ubound( positions, 1 )
      istart2 = Lbound( positions, 2 )
      iend2 = Ubound( positions, 2 )

      if( .not. Grid_has_positions( this ) ) allocate( this%positions(1:3, 1:this%ngpt) )
      this%positions(1:3, 1:this%ngpt) = positions(istart1:iend1, istart2:iend2)
   end subroutine Grid_set_positions


   !********************************************************************
   ! Ions functions
   !********************************************************************


   !! Return .true. if the ions have charges. 
   !! Even though ions should always be initialized with charges, it is not impossible to not have charges. 
   pure function Ions_has_charges( this ) result( hascharges )
      !! Input: ions
      type(ions_t), intent(in) :: this
      !! Output: whether the ions have charges
      logical                  :: hascharges
      !!

      hascharges = Allocated( this%charges )

   end function Ions_has_charges


   !! Return .true. if the ions have positions. 
   !! Even though ions should always be initialized with positions, it is not impossible to not have positions. 
   pure function Ions_has_positions( this ) result( haspositions )
      !! Input: ions
      type(ions_t), intent(in) :: this
      !! Output: whether the ions have positions
      logical                  :: haspositions
      !!

      haspositions = Allocated( this%positions )

   end function Ions_has_positions


   !! Set the charges of ions. Initialize charges if none exist.
   subroutine Ions_set_charges( this, charges )
      !! Input: ions
      type(ions_t), intent(inout) :: this
      !! Input: charges, 1D array of size nions
      real(kind=DP), intent(in)   :: charges(:)
      !! Internal: start and end of input charges
      integer                     :: istart, iend
      !!

      if( Size( charges ) /= this%nions ) &
         & call Error("Ions_set_charges: Size(charges) /= this%nions", Size( charges ), this%nions )

      istart = Lbound( charges, 1 )
      iend = Ubound( charges, 1 )

      if( .not. Ions_has_charges( this ) ) allocate( this%charges(1:this%nions) )
      this%charges(1:this%nions) = charges(istart:iend)

   end subroutine Ions_set_charges


   !! Set the positions of ions. Initialize positions if none exist.
   subroutine Ions_set_positions( this, positions )
      !! Input: ions
      type(ions_t), intent(inout) :: this
      !! Input: charges, 2D array of size (3, nions)
      real(kind=DP), intent(in)   :: positions(:, :)
      !! Internal :: start and end of input positions for each dimension
      integer                     :: istart1, iend1, istart2, iend2
      !!

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
