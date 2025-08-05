      subroutine SetCoup  (nr,ls)
c --- primary function: adjust DOF vector for solver for coupled nodes
c
c  input arguments:
c     nr       (int,sc,in)         - matrix and lv size
c
c  output arguments:
c     ls       (int,ar(nr),inout)  - modified element dof vector
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
c
      integer nr,ls(nr), i

c --- adjust element DOF vector for coupling
      do i = 1,nr
         if (ls(i).eq.0) cycle
         if (btest(M_DofBits(ls(i)),W_CpBit)) then
             ls(i) = M_MasterCp(ls(i))
         endif
      enddo

      return
      end
