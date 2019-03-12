!Define consistent precisions, should be IEE-754 standard
module Precision_mod
   implicit none

   public :: SP, DP, QP

   integer, parameter :: SP = Selected_real_kind( p = 6, r = 37 )
   integer, parameter :: DP = Selected_real_kind( p = 15, r = 307 )
   integer, parameter :: QP = Selected_real_kind( p = 33, r = 4931 )

end module Precision_mod
