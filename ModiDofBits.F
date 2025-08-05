      subroutine ModiDofBits
c --- primary function: modify M_DofBits(D_ULen) vector for boundary conditions
c
#include "impcom.inc"
#include "constants.inc"
#include "dof_info.inc"
c
      integer  foriqr
      integer  n

#if defined(DEBUG)
      call SubBeg('ModiDofBits')
#endif

c --- flag the specified disp in M_DofBits vector
      call AdjForBCs ! assume that constraints are always set

c --- flag the specified forces in M_DofBits vector
      n = foriqr(0,DB_MAXDEFINED)
      if (n .gt. 0)  call AdjForLoad

c --- adjust the M_DofBits vector for couplings
      if (D_NumCp .gt. 0) call AdjForCoup

#if defined(DEBUG)
      call SubEnd('ModiDofBits')
#endif

      return
      end
