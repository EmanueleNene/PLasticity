      subroutine m_swapRow (a,b,n)
c --- primary function: swap the terms in two vectors (with increments)
c
c  input arguments:
c     a        (dp,ar(na,n),inout) - vector a
c     b        (dp,ar(nb,n),inout) - vector b
c     n        (int,sc,in)         - length (second dimension) of vectors a, b
c
c  output arguments:
c     a        (dp,ar(na,n),inout) - vector a
c     b        (dp,ar(nb,n),inout) - vector b
c
#include "impcom.inc"
c
      integer  n,i,j
      double precision  a(*),b(*),con

      j = 1
      do i = 1,n
         if (a(j).ne.b(j)) then
             con = a(j)
             a(j) = b(j)
             b(j) = con
         endif
         j = j + n
      enddo

      return
      end
