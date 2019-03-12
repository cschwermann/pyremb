module Xcpot_mod
   use Precision_mod
   use Output_mod
   use Type_mod
   use Utils_mod !TODO should contain integrate and gradient
  !use libxc wrappers
  #ifdef LIBXC
   use Xcpot_libxc_mod
  #endif
  !alternatively use xcfun wrappers
  #ifdef XCFUN
   use Xcpot_xcfun_mod
  #endif
   implicit none

   public :: Xc_energy, Xc_potential

contains

   !TODO subroutine xc_parse

   subroutine Xc_energy( molecule, funcstring, energy )
      type(molecule_t), intent(in) :: molecule
      character(:), intent(in) :: funcstring
      real(kind=DP), intent(out) :: energy
      !functional id, for exchange, correlation and kinetic in this order
      integer :: funcid(1:3)
      !energy density
      real(kind=DP), allocatable :: enerdens(:), temppot(:)
      character(:), allocatable :: functype

      allocate( enerdens(1:molecule%ngpt) )
      allocate( temppot(1:molecule%ngpt) )

      energy(:) = 0.0_DP

      call Xc_parse(funcstring,funcid,functype)
      select case ( functype )
         case( "LDA" )
            call Xc_ener( funcid, molecule%ngpt, molecule%spinpol, temppot, molecule%density )
            enerdens(:) = enerdens(:) + temppot(:)

            call Xc_ener( funcid, molecule%ngpt, molecule%spinpol, temppot, molecule%density )
            enerdens(:) = enerdens(:) + temppot(:)

            call Xc_ener( funcid, molecule%ngpt, molecule%spinpol, temppot, molecule%density )
            enerdens(:) = enerdens(:) + temppot(:)
         case( "GGA" )
            if( .not. Mol_has_gradient( molecule ) call Gradient( molecule )
            call Xc_ener( funcid, molecule%ngpt, molecule%spinpol, temppot, molecule%density, molecule%gradient )
            enerdens(:) = enerdens(:) + temppot(:)

            call Xc_ener( funcid, molecule%ngpt, molecule%spinpol, temppot, molecule%density, molecule%gradient )
            enerdens(:) = enerdens(:) + temppot(:)

            call Xc_ener( funcid, molecule%ngpt, molecule%spinpol, temppot, molecule%density, molecule%gradient )
            enerdens(:) = enerdens(:) + temppot(:)
         case default
            call Error( "Xc_energy: only LDA and GGA possible!" )
      end select

      deallocate( temppot )

      energy = Integrate( enerdens, molecule%grid )

      deallocate( enerdens )
   end subroutine Xc_energy

   !exchange-correlation-kinetic potential for a molecule and a specified functional
   subroutine Xc_potential( molecule, funcstring, potential )
      type(molecule_t), intent(in) :: molecule
      character(:), intent(in) :: funcstring
      real(kind=DP), intent(out) :: potential(:)
      !functional id, for exchange, correlation and kinetic in this order
      integer :: funcid(1:3)
      real(kind=DP), allocatable :: temppot(:)
      character(:), allocatable :: functype

      potential(:) = 0.0_DP
      if( molecule%spinpol ) then
         allocate( temppot(2 * molecule%ngpt) )
      else
         allocate( temppot(molecule%ngpt )
      end if

      call Xc_parse( funcstring, funcid, functype )
      select case ( functype )
         case( "LDA" )
            call Xc_pot( funcid(1), molecule%ngpt, molecule%spinpol, temppot, molecule%density )
            potential(:) = potential(:) + temppot(:)

            call xc_pot( funcid(2), molecule%ngpt, molecule%spinpol, temppot, molecule%density )
            potential(:) = potential(:) + temppot(:)

            call xc_pot( funcid(3), molecule%ngpt, molecule%spinpol, temppot, molecule%density )
            potential(:) = potential(:) + temppot(:)
         case( "GGA" )
            if( .not. Mol_has_gradient( molecule ) ) call Gradient( molecule )
            call Xc_pot( funcid(1), molecule%ngpt, molecule%spinpol, temppot, molecule%density, molecule%gradient )
            potential(:) = potential(:) + temppot(:)

            call Xc_pot( funcid(2), molecule%ngpt, molecule%spinpol, temppot, molecule%density, molecule%gradient )
            potential(:) = potential(:) + temppot(:)

            call Xc_pot( funcid(3), molecule%ngpt, molecule%spinpol, temppot, molecule%density, molecule%gradient )
            potential(:) = potential(:) + temppot(:)
         case default
            call Error( "Xc_potential: only LDA and GGA possible!" )
      end select

      deallocate( temppot )
   end subroutine Xc_potential

end module Xcpot_mod

!libxc wrapper
#ifdef LIBXC

module Xcpot_libxc_mod
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
      real(kind=DP), intent(out) :: energy(:)
      !internal: de/dsigma in addition to de/drho
      real(kind=DP), allocatable :: vsigma(:)

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
            call Xc_f90_lda_exc( xcfunc, ngpt, density, enerdens )
         !gga, needs gradient, returns also de/dsigma, where sigma is gradrho*gradrho
         !needs some formula to yield the correct potential
         case( XC_FAMILY_GGA )
            if( .not. Present( gradient ) ) &
               call Error( "Xc_pot: GGA functional specified but no gradient given!" )
            call Xc_f90_gga_exc( xcfunc, ngpt, density, gradient, enerdens )
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
            call Xc_f90_lda_vxc( xcfunc, ngpt, density, potential )
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
            call Xc_f90_gga_vxc( xcfunc, ngpt, density, gradient, potential, vsigma )
            !TODO tatsächliche formel einsetzen?
            !entälht vermutlich NOCH MEHR GRADIENTEN
            deallocate( vsigma )
         !all other cases are impossible
         case default
            call Error( "xc_pot: only LDA and GGA possible!" )
      end select

      call Xc_f90_func_end( xcfunc )
   end subroutine Xc_pot

#endif   

end module Xc_pot_libxc_mod
