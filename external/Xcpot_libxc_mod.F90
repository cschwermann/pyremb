module Xcpot_libxc_mod
   use Precision_mod
   use Output_mod
   use Xc_f90_types_m
   use Xc_f90_lib_m
   implicit none

   private :: xcfunc, xcinfo, vmajor, vminor
   public  :: Xc_ener, Xc_pot

   type(xc_f90_pointer_t) :: xcfunc
   type(xc_f90_pointer_t) :: xcinfo
   integer :: vmajor, vminor

contains

   subroutine Xc_ener( functional, ngpt, spinpol, enerdens, density, gradient )
      !functional id, #grid points
      integer, intent(in) :: functional,ngpt
      logical, intent(in) :: spinpol
      real(kind=DP), intent(in) :: density(:)
      real(kind=DP), optional, intent(in) :: gradient(:)
      real(kind=DP), intent(out) :: enerdens(:)

      enerdens(:) = 0.0_DP

      !initialize functional
      if(spinpol) then
         call Xc_f90_func_init( xcfunc, xcinfo, functional, XC_POLARIZED )
      else
         call Xc_f90_func_init( xcfunc, xcinfo, functional, XC_UNPOLARIZED )
      endif

      select case ( Xc_f90_info_family( xcinfo ) )
         !lda, simple
         case( XC_FAMILY_LDA )
            call Xc_f90_lda_exc( xcfunc, ngpt, density(1), enerdens(1) )
         !gga, needs gradient, returns also de/dsigma, where sigma is gradrho*gradrho
         !needs some formula to yield the correct potential
         case( XC_FAMILY_GGA )
            if( .not. Present( gradient ) ) &
               call Error( "Xc_pot: GGA functional specified but no gradient given!" )
            call Xc_f90_gga_exc( xcfunc, ngpt, density(1), gradient(1), enerdens(1) )
         !all other cases are impossible
         case default
            call Error( "Xc_pot: only LDA and GGA possible!" )
      end select

      call Xc_f90_func_end( xcfunc )
   end subroutine Xc_ener

   subroutine Xc_pot( functional, ngpt, spinpol, potential, density, gradient )
      !functional id, #grid points
      integer, intent(in) :: functional,ngpt
      logical, intent(in) :: spinpol
      real(kind=DP), intent(in) :: density(:)
      real(kind=DP), optional, intent(in) :: gradient(:)
      real(kind=DP), intent(out) :: potential(:)
      !internal: de/dsigma in addition to de/drho
      real(kind=DP), allocatable :: vsigma(:)

      potential(:) = 0.0_DP

      !initialize functional
      if( spinpol ) then
         call Xc_f90_func_init( xcfunc, xcinfo, functional, XC_POLARIZED )
      else
         call Xc_f90_func_init( xcfunc, xcinfo, functional, XC_UNPOLARIZED )
      endif

      select case ( Xc_f90_info_family( xcinfo ) )
         !lda, simple
         case( XC_FAMILY_LDA )
            call Xc_f90_lda_vxc( xcfunc, ngpt, density(1), potential(1) )
         !gga, needs gradient, returns also de/dsigma, where sigma is gradrho*gradrho
         !needs some formula to yield the correct potential
         case( XC_FAMILY_GGA )
            if( .not. Present( gradient ) ) &
               call Error( "Xc_pot: GGA functional specified but no gradient given!" )
            if( spinpol ) then
               allocate( vsigma(1:3*ngpt) )
            else
               allocate( vsigma(1:ngpt) )
            endif
            call Xc_f90_gga_vxc( xcfunc, ngpt, density(1), gradient(1), potential(1), vsigma(1) )
            !TODO tatsächliche formel einsetzen?
            !entälht vermutlich NOCH MEHR GRADIENTEN
            deallocate( vsigma )
         !all other cases are impossible
         case default
            call Error( "Xc_pot: only LDA and GGA possible!" )
      end select

      call Xc_f90_func_end( xcfunc )
   end subroutine Xc_pot

   subroutine Xc_parse( funcstring, funcid, functype )
      character(len=*), intent(in) :: funcstring
      integer, intent(out) :: funcid
      character(len=:), allocatable, intent(out) :: functype

      select case( Trim( funcstring ) )
         case( "LDA" )
            funcid = 1
            functype = "LDA"
         case( "PBE" )
            funcid = 2
            functype = "GGA"
         !TODO weitere
         !TODO aufschlüsseln, X, C, K
         case default
            call Error( "Xc_parse: Unknown functional!" )
      end select

   end subroutine Xc_parse

end module Xcpot_libxc_mod
