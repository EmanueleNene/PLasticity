      function vi_last (iv,n)
c --- primary function: find position of last non-zero term in a integer vector
c
c  input arguments:
c     iv       (int,ar(n),in)    - input integer vector
c     n        (int,sc,in)       - size of vector "iv"
c
c  output arguments:
c     vi_last   (int,sc,in)       - position of last non-zero term
c
#include "impcom.inc"
c
      integer  n, iv(n), i, vi_last
c

      do 10 i = n,1,-1
         if (iv(i) .ne. 0)   then
             vi_last = i
             goto 999
         endif
 10   continue

      vi_last = 0
      
 999  return
      end
