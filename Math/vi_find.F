      function  vi_find (iv,ival,n)
c --- primary function: find the location of an integer in a vector
c
#include "impcom.inc"
c
      integer  n,iv(*),ival,vi_find, i
c

      do 10 i = 1,n
         if (iv(i) .eq. ival) then
             vi_find = i
             goto 999
         endif
  10  continue

      vi_find = 0

  999 return
      end
