      subroutine m_diag (a,na,n, eorder)
c --- primary function: diagonalize a matrix
c
c  input arguments:
c       na     (int,sc,in)           - row dimension of matrix a
c       n      (int,sc,in)           - portion to diagonalize
c       eorder (int,sc,in)           - element order
c                                      eorder = 0  -> 2-node elements
c                                      eorder = 1  -> middle node elements
c                                                     (apportion total mass according to
c                                                      fraction of total diagonals)
c
c  output arguments:
c       a      (dp,ar(na,na),inout)  - matrix to diagonalize
c
#include "impcom.inc"
c
      double precision  v_sumi,v_sum

      integer  na,n,i, eorder
      double precision  a(na,na),ccon,dcon
 

      if (eorder.ne.0) then

          dcon = v_sumi (a(1,1),na+1,n)
          if (dcon.gt.0.0d0) then
              ccon = 0.0d0
              do 100 i = 1,n
                 ccon = ccon + v_sum (a(1,i),n)
 100          continue
              ccon = ccon/dcon
              do 200 i = 1,n
                 dcon = ccon*a(i,i)
                 call v_zero (a(1,i),n)
                 a(i,i) = dcon
 200          continue
          endif

      else

          do 300 i = 1,n
             ccon = v_sum (a(1,i),n)
             call v_zero (a(1,i),n)
             a(i,i) = ccon
 300      continue

      endif

      return
      end
