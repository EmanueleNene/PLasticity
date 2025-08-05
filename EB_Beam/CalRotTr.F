      subroutine CalRotTr (rot,tr)
c --- primary function: compute transformation matrix from rotation vector
c
c     ref 1: eqn. 20 (transposed),rankin and brogan, jpvt,108(1986)165-174
c     ref 2: eqn. b.4, simo and vu-quoc, CMANE, 58 (1986), 79-116
c                     (removes singularities at pi and 2*pi)
c
c  input arguments:
c     rot      (dp,ar(3),in)     - rotation parameter in radians
c
c  output arguments:
c     tr       (dp,ar(3,3),out)  - transformation matrix corresponding to rot
c
#include "impcom.inc"
c
      integer i,j
      double precision rot(3),tr(3,3),con,IM(3,3)

c --- 3x3 identity matrix
      data  ((IM(i,j),i=1,3),j=1,3) /
     &  1.0d0,   0.0d0,   0.0d0,
     &  0.0d0,   1.0d0,   0.0d0,
     &  0.0d0,   0.0d0,   1.0d0/

      call BasRot (IM(1,1),tr(1,1),rot(1),3)

c --- transpose
      con = tr(1,2)
      tr(1,2) = tr(2,1)
      tr(2,1) = con
      con = tr(1,3)
      tr(1,3) = tr(3,1)
      tr(3,1) = con
      con = tr(2,3)
      tr(2,3) = tr(3,2)
      tr(3,2) = con

      return
      end
