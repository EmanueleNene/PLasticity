      subroutine MatOrtho (prop, eEl, nTens, dsdeEl, sigma, eTh, siglv)
c --- primary function: compute stress, elastic stiffness matrix
c                       for orthotropic linear elasticity material
c
c  input arguments:
c       prop      (dp,ar(9),in)            - E1,E2,E3,etc. (function of moist)
c       eEl       (dp,ar(nTens),in)        - strainEl_t + strain_Inc
c       nTens     (int,sc,in)              - nDirect + nShear
c       iniData   (dp,ar(*),in)            - initial stress data
c
c
c  output arguments:
c       dsdeEl    (dp,ar(nTens,nTens),out) - elastic stiffness matrix
c       sigma     (dp,ar(nTens),out)       - stress based on elastic predictor
c
c  Mind the order: xx, yy, zz, xy, yz, zx
c
#include "impcom.inc"
c
      integer          nTens
      double precision prop(*), eEl(nTens), dsdeEl(nTens,nTens),
     &                 sigma(nTens), eTh(nTens), siglv(nTens)

      dsdeEl(1,1) = prop(1)
      sigma(1) = prop(1) * eEl(1)
      siglv(1) = prop(1) * eTh(1)

c --- transverse shear
      dsdeEl(2,2) = prop(6)
      dsdeEl(3,3) = prop(4)
      sigma(2) = prop(6) * eEl(2)
      sigma(3) = prop(4) * eEl(3)

      return
      end
