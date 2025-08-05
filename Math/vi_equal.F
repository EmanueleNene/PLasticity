      function vi_equal (v1,v2,n)
c --- primary function: check if two vectors are equal
c
c  input arguments:
c     v1       (int,ar(n),in)    - input integer vector 1
c     v2       (int,ar(n),in)    - input integer vector 2
c     n        (int,sc,in)       - size of vectors
c
c  output arguments:
c     vi_equal  (int,sc,in)      - true if equal
c
#include "impcom.inc"
c
      integer   n, v1(n), v2(n), i
      logical   vi_equal

      do 10 i = 1,n
         if (v1(i).ne.v2(i)) then
             vi_equal = .false.
             goto 999
         endif
 10   continue

      vi_equal = .true.       

 999  return
      end
