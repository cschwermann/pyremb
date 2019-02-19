module xcpot_mod
   use output_mod
   use type_mod
   use utils_mod !TODO should contain integrate and gradient
  !use libxc wrappers
  #ifdef LIBXC
   use xcpot_libxc_mod
  #endif
  !alternatively use xcfun wrappers
  #ifdef XCFUN
   use xcpot_xcfun_mod
  #endif
   implicit none
   public

  contains

   !TODO subroutine xc_parse

   subroutine xc_energy(molecule,funcstring,energy)
      type(T_Molecule),intent(in) :: molecule
      character(:),intent(in) :: funcstring
      real*8,intent(out) :: energy
      !functional id, for exchange, correlation and kinetic in this order
      integer :: funcid(3)
      !energy density
      real*8, allocatable :: enerdens(:),temppot(:)
      character(:),allocatable :: functype

      allocate(enerdens(molecule%ngpt))
      allocate(temppot(molecule%ngpt))

      energy=0.0d0

      call xc_parse(funcstring,funcid,functype)
      select case (functype)
         case("LDA")
            call xc_ener(funcid,molecule%ngpt,molecule%spinpol,temppot,molecule%density)
            enerdens=enerdens+temppot

            call xc_ener(funcid,molecule%ngpt,molecule%spinpol,temppot,molecule%density)
            enerdens=enerdens+temppot

            call xc_ener(funcid,molecule%ngpt,molecule%spinpol,temppot,molecule%density)
            enerdens=enerdens+temppot
         case("GGA")
            if(.not.Mol_has_gradient(molecule) call gradient(molecule)
            call xc_ener(funcid,molecule%ngpt,molecule%spinpol,temppot,molecule%density,molecule%gradient)
            enerdens=enerdens+temppot

            call xc_ener(funcid,molecule%ngpt,molecule%spinpol,temppot,molecule%density,molecule%gradient)
            enerdens=enerdens+temppot

            call xc_ener(funcid,molecule%ngpt,molecule%spinpol,temppot,molecule%density,molecule%gradient)
            enerdens=enerdens+temppot
         case default
            call error("xc_energy: only LDA and GGA possible!")
      end select

      deallocate(temppot)

      energy=integrate(enerdens,molecule%grid)

      deallocate(enerdens)
   end subroutine

   !exchange-correlation-kinetic potential for a molecule and a specified functional
   subroutine xc_potential(molecule,funcstring,potential)
      type(T_Molecule),intent(in) :: molecule
      character(:),intent(in) :: funcstring
      real*8,intent(out) :: potential(:)
      !functional id, for exchange, correlation and kinetic in this order
      integer :: funcid(3)
      real*8,allocatable :: temppot(:)
      character(:),allocatable :: functype

      potential=0.0d0
      if(molecule%spinpol) then
         allocate(temppot(2*molecule%ngpt)
      else
         allocate(temppot(molecule%ngpt)
      endif

      call xc_parse(funcstring,funcid,functype)
      select case (functype)
         case("LDA")
            call xc_pot(funcid(1),molecule%ngpt,molecule%spinpol,temppot,molecule%density)
            potential=potential+temppot

            call xc_pot(funcid(2),molecule%ngpt,molecule%spinpol,temppot,molecule%density)
            potential=potential+temppot

            call xc_pot(funcid(3),molecule%ngpt,molecule%spinpol,temppot,molecule%density)
            potential=potential+temppot
         case("GGA")
            if(.not.Mol_has_gradient(molecule) call gradient(molecule)
            call xc_pot(funcid(1),molecule%ngpt,molecule%spinpol,temppot,molecule%density,molecule%gradient)
            potential=potential+temppot

            call xc_pot(funcid(2),molecule%ngpt,molecule%spinpol,temppot,molecule%density,molecule%gradient)
            potential=potential+temppot

            call xc_pot(funcid(3),molecule%ngpt,molecule%spinpol,temppot,molecule%density,molecule%gradient)
            potential=potential+temppot
         case default
            call error("xc_potential: only LDA and GGA possible!")
      end select

      deallocate(temppot)
   end subroutine

end module

!libxc wrapper
#ifdef LIBXC
module xcpot_libxc_mod
   use xc_f90_types_m
   use xc_f90_lib_m
   implicit none
   private
   type(xc_f90_pointer_t) :: xcfunc
   type(xc_f90_pointer_t) :: xcinfo
   integer :: vmajor, vminor

   public xc_ener
   public xc_pot

  contains

   subroutine xc_ener(functional,ngpt,spinpol,enerdens,density,gradient)
      !functional id, #grid points
      integer,intent(in) :: functional,ngpt
      logical,intent(in) :: spinpol
      real*8,intent(in) :: density(:)
      real*8,optional,intent(in) :: gradient(:)
      real*8,intent(out) :: energy(:)
      !internal: de/dsigma in addition to de/drho
      real*8,allocatable :: vsigma(:)

      enerdens=0.0d0

      !initialize functional
      if(spinpol) then
         call xc_f90_func_init(xcfunc,xcinfo,functional,xc_polarized)
      else
         call xc_f90_func_init(xcfunc,xcinfo,functional,xc_unpolarized)
      endif

      select case (xc_f90_info_family(xcinfo))
         !lda, simple
         case(xc_family_lda)
            call xc_f90_lda_exc(xcfunc,ngpt,density,enerdens)
         !gga, needs gradient, returns also de/dsigma, where sigma is gradrho*gradrho
         !needs some formula to yield the correct potential
         case(xc_family_gga)
            if(.not.present(gradient) &
               call error("xc_pot: GGA functional specified but no gradient given!")
            call xc_f90_gga_exc(xcfunc,ngpt,density,gradient,enerdens)
         !all other cases are impossible
         case default
            call error("xc_pot: only LDA and GGA possible!")
      end select

      call xc_f90_func_end(xcfunc)
   end subroutine

   subroutine xc_pot(functional,ngpt,spinpol,potential,density,gradient)
      !functional id, #grid points
      integer,intent(in) :: functional,ngpt
      logical,intent(in) :: spinpol
      real*8,intent(in) :: density(:)
      real*8,optional,intent(in) :: gradient(:)
      real*8,intent(out) :: potential(:)
      !internal: de/dsigma in addition to de/drho
      real*8,allocatable :: vsigma(:)

      potential=0.0d0

      !initialize functional
      if(spinpol) then
         call xc_f90_func_init(xcfunc,xcinfo,functional,xc_polarized)
      else
         call xc_f90_func_init(xcfunc,xcinfo,functional,xc_unpolarized)
      endif

      select case (xc_f90_info_family(xcinfo))
         !lda, simple
         case(xc_family_lda)
            call xc_f90_lda_vxc(xcfunc,ngpt,density,potential)
         !gga, needs gradient, returns also de/dsigma, where sigma is gradrho*gradrho
         !needs some formula to yield the correct potential
         case(xc_family_gga)
            if(.not.present(gradient) &
               call error("xc_pot: GGA functional specified but no gradient given!")
            if(spinpol) then
               allocate(vsigma(3*ngpt))
            else
               allocate(vsigma(ngpt))
            endif
            call xc_f90_gga_vxc(xcfunc,ngpt,density,gradient,potential,vsigma)
            !TODO tatsächliche formel einsetzen?
            !entälht vermutlich NOCH MEHR GRADIENTEN
            deallocate(vsigma)
         !all other cases are impossible
         case default
            call error("xc_pot: only LDA and GGA possible!")
      end select

      call xc_f90_func_end(xcfunc)
   end subroutine
#endif   
end module
