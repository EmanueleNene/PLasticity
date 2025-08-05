      subroutine m_swapCol (a,b,n)
c --- primary function: swap the terms in two vectors
c
c  input arguments:
c     a        (dp,ar(na,n),inout) - vector a
c     b        (dp,ar(nb,n),inout) - vector b
c     n        (int,sc,in)         - length of vectors a, b
c
c  output arguments:
c     a        (dp,ar(na,n),inout) - vector a
c     b        (dp,ar(nb,n),inout) - vector b
c
#include "impcom.inc"
c
      integer  n,i
      double precision  a(*),b(*),con

      do i = 1,n
         if (a(i) .ne. b(i)) then
             con = a(i)
             a(i) = b(i)
             b(i) = con
         endif
      enddo

      return
      end
