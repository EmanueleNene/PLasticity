      function v_last (v,n)
c --- primary function: find position of last non-zero term in a dp vector
c
c  input arguments:
c     v        (dp,ar(n),in)     - input integer vector
c     n        (int,sc,in)       - size of vector "iv"
c
c  output arguments:
c     v_last    (int,sc,in)      - position of last non-zero term
c
#include "impcom.inc"
c
      integer  n,i,v_last
      double precision  v(*)
c

      do 10  i = n, 1, -1
         if (v(i) .ne. 0.0d0)   then
             v_last = i
             goto 999
         endif
 10   continue

      v_last = 0
      
 999  return
      end
