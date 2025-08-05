      subroutine SwapBytes (num)
c --- primary function: swap bytes
c
#include "impcom.inc"
c
      integer num,dummy,work
      byte bnum(4),bwork(4)
      equivalence (bnum,dummy)
      equivalence (bwork,work)
      
      dummy = num
      bwork(3) = bnum(2)
      bwork(2) = bnum(3)
      bwork(1) = bnum(4)
      bwork(4) = bnum(1)
      num = work

      return
      end
