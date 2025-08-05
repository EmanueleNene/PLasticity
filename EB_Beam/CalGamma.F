      subroutine CalGamma(gam,qlc,tr,alen)
c --- primary function: compute the gamma matrix transpose for beam elements
c
c     compute the gamma matrix which will be used to form the projector matrix 'pte'
c     ref: rankin & nour-omid in computer & struct. vol. 30, pp.257-267, 1988, eqn. (56)
c
c  input arguments:
c     tr     (dp,ar(3,3),in)        - ek current elem transform matrix transpose
c     qlc    (dp,ar(3),in)          - reference vector which records the element
c                                     rotations using mid-triad matrix
c     alen   (dp,sc,in)             - current length of the element
c
c  output arguments:
c     gam    (dp,ar(12,3),out)      - gamma matrix used to form the projector
c
#include "impcom.inc"
c
      double precision
     & gam(12,3),qlc(3),tr(3,3),vect(3),eta,con,alen,con1,SMALL
      parameter (SMALL = 1.0e-10)

      call v_zero(gam(1,1),36)
      call Axv1_v2(tr(1,1),qlc(1),vect(1),3,3)

      if (abs(vect(2)).le.SMALL) then
          eta = 0.0d0
      else
          eta = vect(1)/vect(2)
      endif

      con = 1.0d0/alen
      con1 = eta/alen
      gam(2,3) = -con
      gam(3,1) = con1
      gam(3,2) = con
      gam(4,1) = 1.0d0
      gam(5,1) = -eta
      gam(8,3) = con
      gam(9,1) = -con1
      gam(9,2) = -con

      return
      end
