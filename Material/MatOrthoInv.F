      subroutine MatOrthoInv (prop, sigma, nTens, eEl)
c --- primary function: compute strain from known stress
c                       for orthotropic linear elasticity material
c
c  input arguments:
c       prop      (dp,ar(9),in)            - E1,E2,E3,etc. (function of moist)
c       sigma     (dp,ar(nTens),in)        - stress based on elastic predictor
c       nTens     (int,sc,in)              - nDirect + nShear
c
c  output arguments:
c       eEl       (dp,ar(nTens),out)       - strainEl_t + strain_Inc
c
c  Mind the order: xx, yy, zz, xy, yz, zx
c
#include "impcom.inc"
c
      integer          nTens
      double precision prop(*), eEl(nTens), sigma(nTens)

      eEl(1) = sigma(1)/prop(1)
      eEl(2) = sigma(2)/prop(6)
      eEl(3) = sigma(3)/prop(4)
      
      return
      end
