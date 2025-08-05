      function MakeLongInt (small,large)
c --- primary function: convert two short integers into a long integer
c
c  input arguments:
c     small    (int,sc,in)             - least significant part
c     large    (int,sc,in)             - most significant part
c
c  output arguments:
c     MakeLongInt  (integer*8,sc,out)  - 64 bit integer
c
#include "impcom.inc"
c
      integer*8            MakeLongInt, lwork
      integer*4            small, large, iwork(2)


      equivalence  (lwork,iwork(1))


#if defined (IEEE_SYS)
      iwork(2) = small
      iwork(1) = large
#else
      iwork(1) = small
      iwork(2) = large
#endif
#if !defined (LONGINT2)
      iwork(1) = small
      iwork(2) = 0
#endif

      MakeLongInt = lwork

      return
      end
