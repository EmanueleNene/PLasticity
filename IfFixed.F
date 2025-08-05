      function IfFixed (nr,ls)
c --- primary function: determine if a element is partly constrained
c
c  input arguments:
c     nr       (int,sc,in)       - number of DOF on element
c     ls(nr)   (int,ar(nr),in)   - list of element DOFs
c
c  output arguments:
c     IfFixed  (log,sc,out)      - true if element is constrained
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
c
      integer   nr, ls(nr), i, kdf
      logical   IfFixed

      IfFixed = .true.

      do i = 1,nr
         kdf = ls(i)
         if (kdf .gt. 0) then
c --- it is constrained
            if (btest(M_DofBits(kdf),W_DispBit)) goto 999
         endif
      enddo
      IfFixed = .false.
 
 999  return
      end
