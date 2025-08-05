      subroutine CalHtb(hthb,deba,nnod)
c --- primary function: compute h theta overbar matrix used in beam elements
c
c     compute the h theta overbar matrix which will be used to
c     transform internal forces
c     ref: nour-omid & rankin 'finite rotation analysis ...'
c          CMANE, vol. 93, pp. 353-384, 1991, eqns. (23) & (a.14)
c
c  input arguments:
c     deba   (dp,ar(6*nnod),in)        - deformable rotations displacements
c                                        deba acts on local (rotated) orientation of the element
c     nnod   (int,sc,in)               - number of nodes in an element
c
c  output arguments:
c     hthb   (dp,ar(3,3*nnod),out)     - h theta overbar matrix dependent of theta overbar
c
#include "impcom.inc"
c
      integer nnod,i,j,j1,k,k1,num
      double precision v_dot,small,eta,con,con1,con2,
     &                 hthb(3,3*nnod),deba(6*nnod),wk1(3,3)
c
      num = 9*nnod
      call v_zero (hthb(1,1),num)
      small = 1.0e-8
      j = 3
      do 10 i = 1,nnod
         j1 = (i - 1)*3
         con = v_dot (deba(j+1),deba(j+1),3)
         con1 = sqrt (con)
         if (con1.le.small) then
             do 14 k = 1,3
                hthb(k,k+j1) = 1.0d0
  14         continue 
         else
             con2 = con1*0.50d0
             eta = (sin(con2) - con2*cos(con2))/sin(con2)
             eta = eta/con
             call skewA3x3 (deba(j+1),wk1(1,1))
             do 11 k = 1,3
                do 12 k1 = 1,3
                   wk1(k,k1) = -0.50d0*wk1(k,k1) + eta*deba(k+j)
     &                         *deba(k1+j)
                   if (k.eq.k1)
     &                 wk1(k,k1) = wk1(k,k1) + (1.0d0 - con*eta)
                   hthb(k,k1+j1) = wk1(k,k1)
  12            continue
  11         continue
         endif
  10  j = j+6
c
      return
      end
