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

   !Interpolate the function ref_func, which exists on grid with positions ref_pos
   !onto the grid pos, which results in function values func
   !Imporant: Might need to rescale everything if the integral has to be the same!
   subroutine Shepard_interpolate( ref_pos, ref_func, pos, func )
      real(kind=DP), intent(in) :: ref_pos(:, :), pos(:, :)
      real(kind=DP), intent(in) :: ref_func(:)
      real(kind=DP), intent(out) :: func(:)
      !internal
      integer :: ref_istart, ref_iend, istart, iend, i, j
      real(kind=DP) :: sum_func, sum_weights, value, weight

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

      !Actual interpolation
      !loop over result grid
      RESULTS: do concurrent ( i = istart:iend )
         sum_func = 0.0_DP
         sum_weights = 0.0_DP
         !loop over reference grid
         !TODO OpenMP stuff
         REFERENCES: do concurrent ( j = ref_istart:ref_iend )
            value = ref_func(j)
            weight = Sum( ( pos(:, i) - ref_pos(:, j) ) ** 2 )
            weight = 1.0_DP / weight
            sum_func = sum_func + weight * value
            sum_weights = sum_weights + weight
         end do REFERENCES
         !new value is quotient of sums
         func(i) = sum_func / sum_weights

         !TODO IEEE error handling? important maths here!

      end do RESULTS

   end subroutine Shepard_interpolate
         

   subroutine Gradient( molecule )
      type(molecule_t) :: molecule

   end subroutine Gradient

end module Utils_mod
