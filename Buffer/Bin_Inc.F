      subroutine Bin_Inc (nblk,LongLocL,nLeng,kbfint)
c --- primary function: increment file pointer
c
c  input arguments:
c     nblk     (int,sc,in)       - increment number
c     nLeng    (int,sc,in)       - number of integer*4 words
c     kbfint   (int,sc,in)       - type of data (0=dp, 1=int)
c
c  output arguments:
c     LongLocL (int*8,sc,inout)  - file pointer
c
#include "impcom.inc"
c
      integer*8  LongLocL
      integer  nblk, nLeng, kbfint


c --- increment LongLocL for initial two words
      LongLocL = LongLocL + 2

c --- increment LongLocL for data
#if defined (FULL_64)
      if (kbfint .eq. 0) then
          LongLocL = LongLocL + 2*nLeng
      else
          LongLocL = LongLocL + nLeng
      endif
#else
      LongLocL = LongLocL + nLeng
#endif

c --- increment LongLocL for trailing word
      LongLocL = LongLocL + 1


      return
      end
