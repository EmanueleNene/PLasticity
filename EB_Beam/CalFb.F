      subroutine CalFb(fba,fe,nnod,nr)
c --- primary function: compute f overbar matrix for beam element
c
c      f overbar matrix is used to form the stress stiffening matrix
c      ref: nour-omid&rankin 'finite rotation analysis ...'
c           CMANE, vol,93, pp. 353-384;1991; eqns. (78)
c
c  input arguments:
c     nnod   (int,sc,in)           - number of nodes
c     fe     (dp,ar(6*nnod),in)    - fe force vector (no overbar local and total)
c     nr     (int,sc,in)           - dimension for number of rows in matrix f overbar
c
c  output arguments:
c     fba    (dp,sc,out)           - f overbar matrix used to calculate stress stiffening
c
c
#include "impcom.inc"
c
      integer nnod,i,j,k,k1,num,nr
      double precision
     & fba(nr,3),fe(6*nnod),wk1(3,3)
c
      num = 2*nnod
      j = 0
      do 10 i = 1,num
         call skewA3x3(fe(j+1),wk1(1,1))
         if(mod(i,2).eq.1) then
            do 12 k=1,3
               do 11 k1=1,3
                  fba(k+j,k1) = -wk1(k,k1)
  11           continue
  12        continue
         else
            do 13 k = 1,3
               do 14 k1 = 1,3
                  fba(k+j,k1) = -0.50d0*wk1(k,k1)
  14           continue
  13        continue
         endif
         j = j + 3
  10  continue
      return
      end
