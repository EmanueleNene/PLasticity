      subroutine MatRotTran (a,tr,nr,ntr)
c --- primary function: perform tr-transposed * a * tr
c
c  input arguments:
c     a        (dp,ar(nr,nr),inout) - matrix to be transformed
c     tr       (dp,ar(ntr,ntr),in)  - transformation matrix
c     nr       (int,sc,in)          - dimensioned size of a matrix
c     ntr      (int,sc,in)          - dimensioned size of tr matrix
c
c  output arguments:
c     a        (dp,ar(nr,nr),inout) - transformed matrix
c
#include "impcom.inc"
c
      integer nr,ntr,i,nn,kr
      double precision  a(nr,nr),tr(ntr,ntr)

c --- check for no rotation
      do i = 1,ntr
         if (tr(i,i).lt.1.0d0) goto 100
      enddo

      goto 200

 100  nn = nr/ntr
      kr = 1
      do i = 1,nn
         call RotNode (tr(1,1),a(1,1),kr,nn,ntr,nr)
         kr = kr + ntr
      enddo 

 200  return
      end
