      subroutine MakeShortInt (long64,small,large)
c primary function: convert a long integer into two short integers
c
c  input arguments:
c     long64   (int*8,sc,in)     - 64 bit integer
c
c  output arguments:
c     small    (int,sc,in)       - least significant part
c     large    (int,sc,in)       - most significant part
c
#include "impcom.inc"
c
      integer*8           long64, lwork
      integer*4           small, large, iwork(2)
      equivalence  (lwork,iwork(1))

      lwork = long64

#if defined (IEEE_SYS)
      small = iwork(2)
      large = iwork(1)
#else
      small = iwork(1)
      large = iwork(2)
#endif
#if !defined (LONGINT2)
      small = iwork(1)
      large = 0
#endif

      return
      end
