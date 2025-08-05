      subroutine AxB_C3x3 (a,b,c, n)
c --- primary function: multiply 3 x n matrix times 3 x 3 matrix
c
c  input arguments:
c     a        (dp,ar(3,3),in)   - matrix a (3,3)
c     b        (dp,ar(3,n),in)   - matrix b (3,n)
c     n        (int,sc,in)       - number of columns in matrix b 
c
c  output arguments:
c     c        (dp,ar(3,n),out)  - product matrix c
c
#include "impcom.inc"
c
      integer  n, j
      double precision   a(3,3),b(3,n),c(3,n)

      do j = 1,n
         c(1,j) = a(1,1)*b(1,j) + a(1,2)*b(2,j) + a(1,3)*b(3,j)
         c(2,j) = a(2,1)*b(1,j) + a(2,2)*b(2,j) + a(2,3)*b(3,j)
         c(3,j) = a(3,1)*b(1,j) + a(3,2)*b(2,j) + a(3,3)*b(3,j)
      enddo

      return
      end
