      function vi_compress (ivect,n)
c --- primary function: compress duplicates from integer vector
c --- secondary function: remove zeros
c
c  input arguments:
c     n           (int,sc,in)        - vector lenght
c
c  output arguments:
c     ivect       (int,ar(n),inout)  - vector
c     vi_compress (int,sc,out)       - compressed length of ivect
c
#include "impcom.inc"
c
      integer  vi_compress, vi_find, n,ivect(n), i,j,k

      vi_compress = 0
      if (n .le. 0) goto 999
c --- sort first
      call vi_sort(ivect(1),n)
      
      i = 1
      k = ivect(i)
      j = 2
  10     if (j .gt. n) goto 900
         if (ivect(j) .ne. k) then
            i = i + 1
            ivect(i) = ivect(j)
            k = ivect(i)
         endif
         j = j + 1
         goto 10

 900  k = i

c --- remove zero if it is there
      j = vi_find (ivect(1),0,k)

      if (j.gt.0) then
          vi_compress = k - 1
          do 20 i = j,k - 1
             ivect(i) = ivect(i+1)
  20      continue
      else
          vi_compress = k
      endif
 
 999  return
      end
