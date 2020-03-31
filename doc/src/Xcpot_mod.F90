module Xcpot_mod
   use Precision_mod
   use Output_mod
   use Types_mod
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

   subroutine Xc_energy( molecule, functional, energy )
      type(molecule_t), intent(in) :: molecule
      character(*), intent(in) :: functional
      real(kind=DP), intent(out) :: energy
      !functional id, for exchange, correlation and kinetic in this order
      integer :: funcid
      !energy density
      real(kind=DP), allocatable :: enerdens(:)
      character(:), allocatable :: functype

      allocate( enerdens(1:molecule%ngpt) )

      enerdens(:) = 0.0_DP

      call Xc_parse( functional, funcid, functype )
         
      select case( functype )
         case( "LDA" )
            call Xc_ener( funcid, molecule%ngpt, molecule%spinpol, enerdens, molecule%density )
         case( "GGA" )
            if( .not. Mol_has_gradient( molecule ) ) call Gradient( molecule )
            call Xc_ener( funcid, molecule%ngpt, molecule%spinpol, enerdens, molecule%density, molecule%gradient )
         case default
            call Error( "Xc_energy: only LDA and GGA possible!" )
      end select

      energy = Integrate( enerdens, molecule%grid )

      deallocate( enerdens )
   end subroutine Xc_energy

   !exchange-correlation-kinetic potential for a molecule and a specified functional
   subroutine Xc_potential( molecule, functional, potential )
      type(molecule_t), intent(in) :: molecule
      character(*), intent(in) :: functional
      real(kind=DP), intent(out) :: potential(:)
      !functional id, for exchange, correlation and kinetic in this order
      integer :: funcid
      real(kind=DP), allocatable :: temppot(:)
      character(:), allocatable :: functype

      potential(:) = 0.0_DP

      call Xc_parse( functional, funcid, functype )
      select case( functype )
         case( "LDA" )
            call Xc_pot( funcid, molecule%ngpt, molecule%spinpol, potential, molecule%density )
         case( "GGA" )
            if( .not. Mol_has_gradient( molecule ) ) call Gradient( molecule )
            call Xc_pot( funcid, molecule%ngpt, molecule%spinpol, potential, molecule%density, molecule%gradient )
         case default
            call Error( "Xc_potential: only LDA and GGA possible!" )
      end select

   end subroutine Xc_potential

end module Xcpot_mod
