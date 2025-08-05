      subroutine v_initi (v,inc,n,const)
c --- primary function: set the components of vector v to const by increments
c
c  input arguments:
c     inc      (int,sc,in)         - increment (first dimension) of vector v
c     n        (int,sc,in)         - number of times to be repeated
c     const    (dp,sc,in)          - constant to set components of vector v to
c
c  output arguments:
c     v        (dp,ar(inc,n),out)  - vector v
c
#include "impcom.inc"
c
      integer  inc,n, i,ia
      double precision  v(*),const
  
      ia = 1
      do i = 1,n
         v(ia) = const
         ia = ia + inc
      enddo

      return
      end
