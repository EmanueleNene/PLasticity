      function v_sumi (v,inc,n)
c --- primary function: sum the components of a vector (with increments)
c
c  input arguments:
c     v        (dp,ar(inc,n),in)   - vector v
c     inc      (int,sc,in)         - increment (first dimension) of vector v
c     n        (int,sc,in)         - length (second dimension) of vector v
c
c  output arguments:
c     v_sumi    (dp,sc,out)        - sum of components of vector v
c
#include "impcom.inc"
c
      integer  inc, n, i,ia
      double precision  v_sumi, v(*), con
 
      con = 0.0d0
      ia = 1
      do i = 1,n
         con = con + v(ia)
         ia = ia + inc
      enddo
      v_sumi = con
c
      return
      end
