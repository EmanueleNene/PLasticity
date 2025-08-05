      subroutine BlockDiagMult (blk,v1,v2,nnod)
c --- primary function: multiplication of a block diagonal matrix by an element local vector
c
c  input arguments:
c     blk    (dp,ar(3,*),in)       - block diagonal matrix; each sub-block is 3x3
c     v1     (dp,ar(6*nnod),in)    - input v vector
c     nnod   (int,sc,in)           - number of nodes in an element
c
c  output arguments:
c     v2     (dp,ar(6*nnod),out)   - output vector  v2  =  [blk] * v1
c
#include "impcom.inc"
c      
      integer nnod,i,j,k,num

      double precision blk(3,*),v1(6*nnod),v2(6*nnod)
c
      num = 6*nnod
      call v_move (v1(1),v2(1),num)
c
      do i = 1,nnod
         j = (i-1)*3 + 1
         k = (i-1)*6 + 4
         call Atxv1_v2 (blk(1,j),v1(k),v2(k),3,3)
      enddo
c      
      return
      end
