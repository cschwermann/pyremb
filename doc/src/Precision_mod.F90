!! Author:  Christian Schwermann
!! E-mail:  c.schwermann@wwu.de
!! Date:    15/03/2019
!! Project: PEREMB 
!! File:    Precision_mod.F90 
!! Copyright: © 2019 Christian Schwermann, ALL RIGHTS RESERVED
!!
!!***********************************************************************
!!
!! Consistent definition of precisions, should be IEE-754 standard
!!
!!***********************************************************************
module Precision_mod
   implicit none

   ! Single  double and quadruple precision
   public :: SP, DP, QP

   !! Single precision
   integer, parameter :: SP = Selected_real_kind( p = 6, r = 37 )
   !! Double precision
   integer, parameter :: DP = Selected_real_kind( p = 15, r = 307 )
   !! Quadruple precision
   integer, parameter :: QP = Selected_real_kind( p = 33, r = 4931 )

end module Precision_mod
