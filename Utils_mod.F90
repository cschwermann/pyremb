!Utility module, contains useful subroutines like integration, etc.
module Utils_mod
   use Precision_mod
   use Types_mod
   implicit none

   public :: Integrate, Shepard_interpolate
   public :: Gradient

contains

   function Integrate( array, grid ) result( integral )
      real(kind=DP), intent(in) :: array(:)
      type(grid_t), optional, intent(in) :: grid
      real(kind=DP) :: integral
      integer :: i, istart, iend

      integral = 0.0_DP

      istart = Lbound( array, 1 )
      iend = Ubound( array, 1 )

      if( Present( grid ) ) then
         if( Size( array ) /= grid%ngpt ) &
            & call Error( "Integrate: Size( array ) /= grid%ngpt", Size( array ), grid%ngpt )

         if( Grid_has_weights( grid ) ) then
            if( istart /= Lbound( grid%weights, 1 ) ) &
               & call Error( "Integrate: Lbound( array ) /= LBound( grid%weights )", istart, LBound( grid%weights, 1 ) )
            if( iend /= Ubound( grid%weights, 1 ) ) &
               & call Error( "Integrate: Ubound( array ) /= UBound( grid%weights )", iend, UBound( grid%weights, 1 ) )
            !TODO OpenMP stuff
            do i = istart, iend
               integral = integral + grid%weights(i) * array(i)
            end do
         else
            do i = istart, iend
               integral = integral + array(i)
            end do
         end if
      else 
         do i = istart, iend
            integral = integral + array(i)
         end do
      end if

   end function Integrate

   !Interpolate all values of a molecule onto some grid
   subroutine Mol_interpolate( molecule, grid )
      type(molecule_t), intent(inout) :: molecule
      type(grid_t), intent(in) :: grid
      integer :: ngpt
      real(kind=DP), allocatable :: newvalues(:)

      !Catch some errors
      if( .not. Mol_has_grid( molecule ) )&
         & call Error( "Mol_interpolate: Molecule has no grid!" )

      !TODO Calculate positions from cell, if cell exists
      if( .not. Grid_has_positions( molecule%grid ) ) &
         & call Error( "Mol_interpolate: Molecule%grid has no positions!" )
      if( .not. Grid_has_positions( grid ) ) &
         & call Error( "Mol_interpolate: Grid has no positions!" )

      !Set new grid points, so all the "Mol_set" later on works
      !New values always have dimension of grid
      ngpt = grid%ngpt
      molecule%ngpt = ngpt
      allocate( newvalues(ngpt) )

      !for spin polarized case, all the arrays are different sizes
      if( molecule%spinpol ) then
         !interpolate density
         if( Mol_has_density( molecule ) ) then
            call Shepard_interpolate( molecule%grid%positions, &
            & molecule%density(1:ngpt), grid%positions, newvalues, .true. )
            call Mol_set_density_a( molecule, newvalues )
            call Shepard_interpolate( molecule%grid%positions, &
            & molecule%density(ngpt+1:2*ngpt), grid%positions, newvalues, .true. )
            call Mol_set_density_b( molecule, newvalues )
         endif
         !interpolate gradient
         if( Mol_has_gradient( molecule ) ) then
            call Shepard_interpolate( molecule%grid%positions, &
            & molecule%gradient(1:ngpt), grid%positions, newvalues, .true. )
            call Mol_set_gradient_aa( molecule, newvalues )
            call Shepard_interpolate( molecule%grid%positions, &
            & molecule%gradient(ngpt+1:2*ngpt), grid%positions, newvalues, .true. )
            call Mol_set_gradient_ab( molecule, newvalues )
            call Shepard_interpolate( molecule%grid%positions, &
            & molecule%gradient(2*ngpt+1:3*ngpt), grid%positions, newvalues, .true. )
            call Mol_set_gradient_bb( molecule, newvalues )
         endif

      else
         !interpolate density
         if( Mol_has_density( molecule ) ) then
            call Shepard_interpolate( molecule%grid%positions, &
            & molecule%density, grid%positions, newvalues, .true. )
            call Mol_set_density( molecule, newvalues )
         endif
         !interpolate gradient
         if( Mol_has_gradient( molecule ) ) then
            call Shepard_interpolate( molecule%grid%positions, &
            & molecule%gradient, grid%positions, newvalues, .true. )
            call Mol_set_gradient( molecule, newvalues )
         endif

      end if

      !interpolate nuclear potential
      if( Mol_has_vnuc( molecule ) ) then
         call Shepard_interpolate( molecule%grid%positions, &
         & molecule%vnuc, grid%positions, newvalues )
         call Mol_set_vnuc( molecule, newvalues )
      endif

      !Finally, set molecular grid to new grid
      call Mol_set_grid( molecule, grid )

   end subroutine Mol_interpolate

   !Interpolate the function ref_func, which exists on grid with positions ref_pos
   !onto the grid pos, which results in function values func
   !Imporant: Might need to rescale everything if the integral has to be the same!
   subroutine Shepard_interpolate( ref_pos, ref_func, pos, func, texp )
      real(kind=DP), intent(in) :: ref_pos(:, :), pos(:, :)
      real(kind=DP), intent(in) :: ref_func(:)
      real(kind=DP), intent(out) :: func(:)
      logical, optional, intent(in) :: texp
      !internal
      integer :: ref_istart, ref_iend, istart, iend, i, j, izero
      real(kind=DP) :: sum_func, sum_weights, value, weight
      real(kind=DP), allocatable :: distances(:)

      !reference start and end index
      ref_istart = Lbound( ref_func, 1 )
      ref_iend = Ubound( ref_func, 1 )
      if( ref_istart /= LBound( ref_pos, 2 ) ) &
         & call Error( "Shepard_Interpolate: LBound( ref_func ) /= LBound( ref_pos )", ref_istart, LBound( ref_pos, 2 ) )
      if( ref_iend /= UBound( ref_pos, 2 ) ) &
         & call Error( "Shepard_Interpolate: UBound( ref_func ) /= UBound( ref_pos )", ref_iend, UBound( ref_pos, 2 ) )

      !result start and end index
      istart = Lbound( func, 1 )
      iend = Ubound( func, 1 )
      if( istart /= LBound( pos, 2 ) ) &
         & call Error( "Shepard_Interpolate: LBound( func ) /= LBound( pos )", istart, LBound( pos, 2 ) )
      if( iend /= UBound( pos, 2 ) ) &
         & call Error( "Shepard_Interpolate: UBound( func ) /= UBound( pos )", ref_iend, UBound( pos, 2 ) )

!TODO add check for first dimension of positions
!      if( istart /= LBound( pos, 2 ) ) &
!         & call Error( "Shepard_Interpolate: LBound( func ) /= LBound( pos )", istart, LBound( pos, 2 ) )

      !allocate distances, they have the length of the reference value array
      allocate( distances(ref_istart:ref_iend) )

      !Actual interpolation
      !loop over result grid
      RESULTS: do concurrent ( i = istart:iend )
         sum_func = 0.0_DP
         sum_weights = 0.0_DP
         izero = -1
         !calculate distances beforehand, as we have to know if one of them is zero
         DISTANCE: do concurrent ( j = ref_istart:ref_iend )
            value = Sum( ( pos(:, i) - ref_pos(:, j) ) ** 2 )
            distances(j) = value
            !this doesnt need to be safe, as the ref_func value shoudl be the same wherever the distance is zero!
            if( value <= 1.0E-6_DP ) izero = j
         end do DISTANCE
         if( Present( texp ) ) then
            if( texp ) distances(:) = Sinh( Sqrt( distances(:) ) )
         end if
         !find minimum to check if the grid point is known
         !It might be worth it to sort here instead of Minloc, that way one could later on
         !sum over N values nearest to reference grid point
         if( izero > 0 ) then 
            !if grid point is known, just use that value
            func(i) = ref_func(izero)
         else
            !loop over reference grid
            !TODO OpenMP stuff
            REFERENCES: do concurrent ( j = ref_istart:ref_iend )
               value = ref_func(j)
               !weight = Sum( ( pos(:, i) - ref_pos(:, j) ) ** 2 )
               weight = 1.0_DP / distances(j)
               sum_func = sum_func + weight * value
               sum_weights = sum_weights + weight
            end do REFERENCES
            !new value is quotient of sums
            func(i) = sum_func / sum_weights
         end if

         !TODO IEEE error handling? important maths here!

      end do RESULTS

   end subroutine Shepard_interpolate
         

   subroutine Gradient( molecule )
      type(molecule_t) :: molecule

   end subroutine Gradient

end module Utils_mod
