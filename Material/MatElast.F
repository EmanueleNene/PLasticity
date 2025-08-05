      subroutine MatElast (nTens, prop, strainEl, stress,
     &                     dsdeEl, strainTh, siglv)
c --- primary function: compute stress, elastic stiffness matrix,
c                       for linear elasticity material
c
c  input arguments:
c      nTens     (int,sc,in)               - nDirect + nShear
c      strainEl  (dp,ar(nTens),in)         - strain_El at t + strain_Inc
c
c  output arguments:
c      dsdeEl     (dp,ar(nTens,nTens),out) - elastic stiffness matrix
c      stress     (dp,ar(nTens),out)       - stress based on el predictor
c      siglv      (dp,ar(*),out)           - hygro stress vector
c
#include "impcom.inc"
c
      integer          nTens

      double precision strainEl(nTens), dsdeEl(nTens,nTens),
     &                 stress(nTens), strainTh(nTens),
     &                 siglv(nTens), prop(21)

      call v_zero (dsdeEl(1,1),nTens*nTens)
      call v_zero (stress(1),nTens)
      call v_zero (siglv(1),nTens)

c --- orthotropic elasticity
      call MatOrtho (prop(1), strainEl(1), nTens,
     &               dsdeEl(1,1), stress(1),
     &               strainTh(1), siglv(1))
c

      return
      end
