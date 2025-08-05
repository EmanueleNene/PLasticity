      subroutine ElasticDriver (nTens, elemId, matId, moist,
     &                          matMW, matMWIndex, cJac, prop,
     &                          tMatProp, khygro, siglv, kMatPoint)
c --- primary function: performs elastic constitutive calculations
c
c  input arguments:
c      elemId    (int,sc,in)        - element number
c      matId     (int,sc,in)        - material number
c      moist     (dp,sc,in)         - moisture
c      matMWIndex(int,ar(*),in)     - array indexing into matMW
c      khygro    (int,sc,in)        - hygro flag
c
c  output arguments:
c      matMW    (dp,ar(*),inout)    - mat record in work array
c      cJac     (dp,ar(nt,*),out)   - consistent Jacobian
c      tMatProp (dp,ar(*),out)      - mat prop for tr shear calc
c      siglv    (dp,ar(4),out)      - hygro load vector
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "solu_info.inc"
c
      double precision HALF, ONE, TOLER
      parameter       (HALF  = 0.5d0,
     &                 ONE   = 1.0d0,
     &                 TOLER = 1.0d-10)
      integer          nTens, elemId, matId, kMatPoint,
     &                 matMWIndex(*), khygro
      double precision moist, matMW(*), cJac(nTens,*),
     &                 tMatProp(*), siglv(*)
      double precision prop(13), g, con, tmpwrk(6)

c --- orthotropic elastic properties
      
      if (kMatPoint.eq.1 .or. khygro.gt.0) then
          call GetElasticProps (elemId, matId, moist, prop(1))

          tMatProp(1) = prop(1)
          tMatProp(2) = prop(5)
          tMatProp(3) = prop(6)
          tMatProp(4) = prop(4)
          tMatProp(5) = prop(13)

          con = abs(prop(1) - prop(2)) + abs(prop(2)- prop(3)) +
     &          abs(prop(4) - prop(5)) + abs(prop(5)- prop(6))
          con = con / (prop(1) + prop(2) + prop(3)) +
     &          abs(prop(7) - prop(8)) + abs(prop(8)- prop(9))
          if (con.lt.TOLER) then
              con = HALF * prop(1) / prop(4) - (ONE + prop(7))
              if (abs(con).lt.TOLER) then
                  g = HALF * prop(1) / (ONE + prop(7))
                  tMatProp(2) = g
                  tMatProp(3) = g
                  tMatProp(4) = g
              endif
          endif
      else
          call v_zero (prop(1),13)
          prop(1) = tMatProp(1)
          prop(5) = tMatProp(2)
          prop(6) = tMatProp(3)
          prop(4) = tMatProp(4)
          prop(13) = tMatProp(5)
      endif

      call MatElast (nTens,prop(1),
     &               matMW(matMWIndex(W_iStrain)),
     &               matMW(matMWIndex(W_iStress)),
     &               matMW(matMWIndex(W_iMatMatrix)),
     &               matMW(matMWIndex(W_iStrainHyg)),
     &               tmpwrk(1))
     
c
      call v_move (matMW(matMWIndex(W_iMatMatrix)),
     &             cJac(1,1), nTens*nTens)
c
      if (khygro.eq.1) then
          if (A_KeyFirstIter.eq.1)  siglv(1) = tmpwrk(1)
      endif

      return
      end
