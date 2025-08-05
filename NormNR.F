      subroutine NormNR (rnorm,rnval,ityp)
c --- primary function: accumulate rnorm by rnval according to ityp.
c                       rnorm is probably the norm of newton-raphson restoring force
c
c  input arguments:
c     rnorm     (dp,sc,inout)      - previously accumulated norm.
c     rnval     (dp,sc,in)         - value to be accumulated
c     ityp      (int,sc,in)        - key to indicate norm used.
c                                    ityp=1: l1-norm or abssolute value norm
c                                    ityp=2: l2-norm or euclidean norm
c
c  output arguments:
c     rnorm     (dp,sc,inout)      - updated norm
c
#include "impcom.inc"
c
      double precision rnorm, rnval
      integer ityp

      if (ityp.eq.2) then
c --- l2 norm
         rnorm = rnorm + rnval*rnval
      elseif (ityp.eq.1) then
c --- l1 norm
         rnorm = rnorm + abs(rnval)
      endif

      return
      end
