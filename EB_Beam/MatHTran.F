      subroutine MatHTran (a,h,nr,ntr)
c --- primary function: perform h-transpose * a * h
c
c  input arguments:
c     a        (dp,ar(nr,nr),inout) - matrix to be transformed
c     h        (dp,ar(ntr,ntr),in)  - h transformation matrix
c     nr       (int,sc,in)          - dimensioned size of a matrix
c     ntr      (int,sc,in)          - dimensioned size of h matrix
c
c  output arguments:
c     a        (dp,ar(nr,nr),inout) - transformed matrix
c
#include "impcom.inc"
c
      integer nr,ntr,i,j,nn,kr
      double precision  a(nr,nr),h(ntr,*)

      nn = nr/ntr
      kr = 4
      j = 1
      do i = 2,nn,2
         call RotNode (h(1,j),a(1,1),kr,nn,ntr,nr)
         kr = kr + ntr*2
         j = j + ntr
      enddo

      return
      end
