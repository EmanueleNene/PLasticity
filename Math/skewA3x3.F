      subroutine skewA3x3(v,t)
c --- primary function: compute skew-symmetric matrix using vector v
c
c  input arguments:
c     v   (dp,ar(3),in)    - axial components of the skew-symmetric matix
c
c  output arguments:
c     t   (dp,ar(3,3),in)  - the skew-symmetric matrix
c
#include "impcom.inc"
c
      double precision v(3),t(3,3)
c
      call v_zeroi (t(1,1),4,3)
      t(1,2) = -v(3)
      t(1,3) =  v(2)
      t(2,1) =  v(3)
      t(2,3) = -v(1)
      t(3,1) = -v(2)
      t(3,2) =  v(1)
      return
      end
