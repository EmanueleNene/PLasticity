      function PowerOf2 (LengthIn,iBitSize)
c --- primary function: find the power of two equal or less than LengthIn
c
c  input arguments:
c
c     LengthIn    (int,sc,in)       - number to process
c     iBitSize    (int,sc,in)       - number of bits for divide
c
#include "impcom.inc"
c
      integer LengthIn, iBitSize

      integer  PowerOf2, i, j

      if (LengthIn.le.0) then
          PowerOf2 = 0
          iBitSize = 0
          goto 999
      endif

      i = 2
      do j = 1,30
         if (LengthIn .lt. i) goto 999
         PowerOf2 = i
         iBitSize = j
         i = i * 2
      enddo

  999 continue

      return
      end
