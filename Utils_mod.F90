!! Author:  Christian Schwermann
!! E-mail:  c.schwermann@wwu.de
!! Date:    15/03/2019
!! Project: PEREMB 
!! File:    Utils_mod.F90 
!! Copyright: © 2019 Christian Schwermann, ALL RIGHTS RESERVED
!!
!!***********************************************************************
!!
!! Utility module, contains generally useful subroutines and functions
!! like integration, interpolation and other "light" math
!!
!!***********************************************************************
module Utils_mod
   use Precision_mod ! IEEE precision
   use Types_mod     ! molecule and grid
   use Output_mod    ! Error and warning
   implicit none

   ! Basic integration, Shepard interpolation and general interpolation for molecules
   public :: Integrate, Shepard_interpolate, Mol_interpolate
   ! Grid volume, cross product, grid positions
   public :: Grid_volume, Cross, Grid_pos_from_cell
   ! Arbitrary gradients
   public :: Gradient


contains


   !!********************************************************************
   !! Integration of arbitrary functions \(a_i\)  on arbitrary grids.
   !!
   !! This is simply a sum, i. e. \( \sum_{i=1}^{N_{gpt}} w_i \cdot a_i \)
   !! when weights \(w_i\) are given and \( \sum_{i=1}^{N_{gpt}} a_i \) else.
   !! 
   !! Here array is \(a_i\), grid%weights is \(w_i\), ngpot is \(N_{gpt}\).
   !!********************************************************************
   function Integrate( array, grid ) result( integral )
      !! Input: any array
      real(kind=DP), intent(in) :: array(:)
      !! Input: grid, has to have the same number of points as the array
      type(grid_t), optional, intent(in) :: grid
      !! Output: the calculated integral
      real(kind=DP) :: integral
      !! Internal: loop index, start and end of array
      integer :: i, istart, iend
      !! Internal: volume of periodic cell
      real(kind=DP) :: volume
      !!
      !!@TODO Specific method for cells @ENDTODO

      integral = 0.0_DP

      istart = Lbound( array, 1 )
      iend = Ubound( array, 1 )

      ! Calculation depends on the grid having weights, a cell or nothing
      if( Present( grid ) ) then
         if( Size( array ) /= grid%ngpt ) &
            & call Error( "Integrate: Size( array ) /= grid%ngpt", Size( array ), grid%ngpt )

         if( Grid_has_weights( grid ) ) then
            if( istart /= Lbound( grid%weights, 1 ) ) &
               & call Error( "Integrate: Lbound( array ) /= LBound( grid%weights )", istart, LBound( grid%weights, 1 ) )
            if( iend /= Ubound( grid%weights, 1 ) ) &
               & call Error( "Integrate: Ubound( array ) /= UBound( grid%weights )", iend, UBound( grid%weights, 1 ) )
            integral = Sum( grid%weights(:) * array(:) )
         else ! no weights
            integral = Sum( array )
            
            ! weight with cell volume
            if( Grid_has_cell( grid ) ) then
               volume = Grid_volume( grid )
               integral = integral / volume
            end if
         end if
      else ! no grid
         integral = Sum( array )
      end if

   end function Integrate


   !!********************************************************************
   !! Calculation of the volume of the cell associated with a given grid.
   !!********************************************************************
   function Grid_volume( grid ) result( volume )
      !! Input: grid, has to have a cell
      type(grid_t), intent(in) :: grid
      !! Output: the calculated volume
      real(kind=DP) :: volume
      !! Internal: cell
      real(kind=DP) :: cell(1:3, 1:3)
      !!

      ! Calculation depends on the grid having weights, a cell or nothing
      if( .not. Grid_has_cell( grid ) ) &
         & call Error( "Grid_volume: Need a cell to calculate volume!" )

      cell(1:3, 1:3) = grid%cell(1:3, 1:3)

      volume = Dot_product( Cross( cell(1, :), cell(2, :) ), cell(3, :) )

   end function Grid_volume


   !!********************************************************************
   !! Calculation of the cross product of two vectors
   !!********************************************************************
   pure function Cross( a, b ) result( c )
      !! Input: two 3D vectors
      real(kind=DP), intent(in) :: a(1:3), b(1:3)
      !! Output: one 3D vector
      real(kind=DP) :: c(1:3)
      !!

      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)

   end function Cross


   !!********************************************************************
   !! Calculation of grid positions from cell vectors and grid points.
   !!********************************************************************
   subroutine Grid_pos_from_cell( grid )
      !! Input: grid, has to have a cell
      type(grid_t), intent(inout) :: grid
      !! Internal: calculated positions
      real(kind=DP), allocatable :: positions(:, :)
      !! Internal: loop indices
      integer :: ix, iy, iz, ii
      !!

      ! Calculation depends on the grid having weights, a cell or nothing
      if( .not. Grid_has_periodicity( grid ) ) &
         & call Error( "Grid_pos_from_cell: Need a cell and points in each direction to calculate positions!" )

      allocate( positions(1:3, 1:grid%ngpt) )

      ii = 0

      !! @NOTE here, x is the fastest dimension,
      !! this is different for cube files, where z is fastest,
      !! thus, reading cubes and converting into our format will be slow @ENDNOTE
      do concurrent ( ix = 1:grid%ngptx )
         do concurrent ( iy = 1:grid%ngpty )
            do concurrent ( iz = 1:grid%ngptz )
               ii = ii + 1
               positions(:, ii) = ix * grid%cell(1, :) + iy * grid%cell(2, :) + iz * grid%cell(3, :)
            end do
         end do
      end do

      call Grid_set_positions( grid, positions )

      deallocate( positions )

   end subroutine Grid_pos_from_cell


   !!********************************************************************
   !! Interpolation of all functions in a molecule onto some grid.
   !!
   !! This subroutine interpolates density, gradient and nuclear potential
   !! of a molecule onto a given grid.
   !! An implementation of Shepards method is used.
   !!********************************************************************
   subroutine Mol_interpolate( molecule, ingrid )
      !! Input: molecule, will be modified
      type(molecule_t), intent(inout) :: molecule
      !! Input: grid
      type(grid_t), intent(in) :: ingrid
      !! Internal: grid, because the input grid shall not be modified
      type(grid_t)             :: grid
      !! Internal: number of grid points of given grid
      integer :: ngpt
      !! Internal: new / interpolated function values, size = ngpt
      real(kind=DP), allocatable :: newvalues(:)
      !! 

      grid = ingrid

      ! Catch some errors
      if( .not. Mol_has_grid( molecule ) )&
         & call Error( "Mol_interpolate: Molecule has no grid!" )

      ! Can calculate positions from cell, if cell exists
      if( .not. Grid_has_positions( molecule%grid ) ) then
         if( Grid_has_periodicity( molecule%grid ) ) then
            call Grid_pos_from_cell( molecule%grid )
         else
            call Error( "Mol_interpolate: Molecule%grid has no positions!" )
         end if
      end if

      if( .not. Grid_has_positions( grid ) ) then
         if( Grid_has_periodicity( grid ) ) then
            call Grid_pos_from_cell( grid )
         else
            call Error( "Mol_interpolate: Grid has no positions!" )
         end if
      end if

      ! Set new grid points, so all the "Mol_set" later on works
      ! New values always have dimension of grid
      ngpt = grid%ngpt
      molecule%ngpt = ngpt
      allocate( newvalues(ngpt) )

      ! For spin polarized case, all the arrays are different sizes
      if( molecule%spinpol ) then
         ! Interpolate density
         if( Mol_has_density( molecule ) ) then
            call Shepard_interpolate( molecule%grid%positions, &
            & molecule%density(1:ngpt), grid%positions, newvalues )
            call Mol_set_density_a( molecule, newvalues )
            call Shepard_interpolate( molecule%grid%positions, &
            & molecule%density(ngpt+1:2*ngpt), grid%positions, newvalues )
            call Mol_set_density_b( molecule, newvalues )
         endif
         ! Interpolate gradient
         if( Mol_has_gradient( molecule ) ) then
            call Shepard_interpolate( molecule%grid%positions, &
            & molecule%gradient(1:ngpt), grid%positions, newvalues )
            call Mol_set_gradient_aa( molecule, newvalues )
            call Shepard_interpolate( molecule%grid%positions, &
            & molecule%gradient(ngpt+1:2*ngpt), grid%positions, newvalues )
            call Mol_set_gradient_ab( molecule, newvalues )
            call Shepard_interpolate( molecule%grid%positions, &
            & molecule%gradient(2*ngpt+1:3*ngpt), grid%positions, newvalues )
            call Mol_set_gradient_bb( molecule, newvalues )
         endif

      else
         ! Interpolate density
         if( Mol_has_density( molecule ) ) then
            call Shepard_interpolate( molecule%grid%positions, &
            & molecule%density, grid%positions, newvalues )
            call Mol_set_density( molecule, newvalues )
         endif
         ! Interpolate gradient
         if( Mol_has_gradient( molecule ) ) then
            call Shepard_interpolate( molecule%grid%positions, &
            & molecule%gradient, grid%positions, newvalues )
            call Mol_set_gradient( molecule, newvalues )
         endif

      end if

      ! Interpolate nuclear potential, doesnt depend on spin
      if( Mol_has_vnuc( molecule ) ) then
         call Shepard_interpolate( molecule%grid%positions, &
         & molecule%vnuc, grid%positions, newvalues )
         call Mol_set_vnuc( molecule, newvalues )
      endif

      ! Finally, set molecular grid to new grid
      call Mol_set_grid( molecule, grid )

   end subroutine Mol_interpolate


   !!********************************************************************
   !! Shepard interpolation of a 3D function \(f_j\) given on positions \(\pmb{r}_{ref,j}\)
   !! onto a grid defined by positions \(\pmb{r}_i\).
   !!
   !! This is the basic, global Shepard method:
   !! $$ Q_i = \frac{ \sum_{j=1}^{M} w_i(\pmb{r}_{ref,j})\cdot f_j }{ \sum_{j=1}^{M} w_i(\pmb{r}_{ref,j}) } $$
   !! with
   !! $$ w_i(\pmb{r}_j) = \frac{1}{d_{i,j}^2} $$
   !! and 
   !! $$ d_{i,j}^2 = |\pmb{r}_i - \pmb{r}_{ref,j}|^2 . $$
   !!
   !! ref_pos is \(\pmb{r}_{ref,j}\), ref_func is \(f_j\), pos is \(\pmb{r}_i\) and func is \(Q_i\).
   !! 
   !! Ref.: Shepard D (1968) A two-dimensional interpolation function for irregularly spaced data 
   !!       Proc. 23rd Nat. Conf. ACM 517-523 Brandon/Systems Press Inc., Princeton 
   !! DOI: [https://dx.doi.org/10.1145/800186.810616](https://dx.doi.org/10.1145/800186.810616)
   !!
   !! @NOTE Might need to rescale everything if the integral has to be the same! @ENDNOTE
   !!********************************************************************
   subroutine Shepard_interpolate( ref_pos, ref_func, pos, func )
      !! Input: reference position and target position
      real(kind=DP), intent(in) :: ref_pos(:, :), pos(:, :)
      !! Input: reference function
      real(kind=DP), intent(in) :: ref_func(:)
      !! Output: interpolated function
      real(kind=DP), intent(out) :: func(:)
      !! Internal: start and end of reference function, start and end of result function, 
      !! loop indices, index of zero distance
      integer :: ref_istart, ref_iend, istart, iend, i, j, izero
      !! Internal: sum w*f, sum w, one weight
      real(kind=DP) :: sum_func, sum_weights, weight
      !! 

      ! Reference start and end index
      ref_istart = Lbound( ref_func, 1 )
      ref_iend = Ubound( ref_func, 1 )
      if( ref_istart /= LBound( ref_pos, 2 ) ) &
         & call Error( "Shepard_Interpolate: LBound( ref_func ) /= LBound( ref_pos )", ref_istart, LBound( ref_pos, 2 ) )
      if( ref_iend /= UBound( ref_pos, 2 ) ) &
         & call Error( "Shepard_Interpolate: UBound( ref_func ) /= UBound( ref_pos )", ref_iend, UBound( ref_pos, 2 ) )

      ! Result start and end index
      istart = Lbound( func, 1 )
      iend = Ubound( func, 1 )
      if( istart /= LBound( pos, 2 ) ) &
         & call Error( "Shepard_Interpolate: LBound( func ) /= LBound( pos )", istart, LBound( pos, 2 ) )
      if( iend /= UBound( pos, 2 ) ) &
         & call Error( "Shepard_Interpolate: UBound( func ) /= UBound( pos )", ref_iend, UBound( pos, 2 ) )

!TODO add check for first dimension of positions
!      if( istart /= LBound( pos, 2 ) ) &
!         & call Error( "Shepard_Interpolate: LBound( func ) /= LBound( pos )", istart, LBound( pos, 2 ) )

      ! Actual interpolation
      ! Loop over result grid
      RESULTS: do concurrent ( i = istart:iend )
         sum_func = 0.0_DP
         sum_weights = 0.0_DP
         ! Loop over reference grid
         !TODO OpenMP stuff
         REFERENCES: do j = ref_istart, ref_iend
            weight = Sum( ( pos(:, i) - ref_pos(:, j) ) ** 2 )
            ! If grid point is known (i.e. distance is 0), just use that value
            ! There is no faster way to do this, i tried
            if ( weight <= 1.0E-6_DP ) then
               sum_func = ref_func(j)
               sum_weights = 1.0_DP
               exit
            end if
            weight = 1.0_DP / weight
            sum_func = sum_func + weight * ref_func(j)
            sum_weights = sum_weights + weight
         end do REFERENCES
         ! New value is quotient of sums
         func(i) = sum_func / sum_weights

         !TODO IEEE error handling? important maths here!

      end do RESULTS

   end subroutine Shepard_interpolate
         


   !!********************************************************************
   !! In future, this will calculate gradients of some function on arbitrary grids
   !!********************************************************************
   subroutine Gradient( molecule )
      ! Input: molecule, will be modifies
      type(molecule_t) :: molecule

   end subroutine Gradient


end module Utils_mod
