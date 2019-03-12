!Utility module, contains useful subroutines like integration, etc.
module Utils_mod
   use Precision_mod
   use Types_mod
   implicit none

   public :: Integrate

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

   subroutine Gradient( molecule )
      type(molecule_t) :: molecule

   end subroutine Gradient

end module Utils_mod
