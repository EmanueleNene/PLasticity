      subroutine FindRoot2 (coef, roots, numRoots)
c --- primary function: solve coef(1)*x^2 + coef(2)*x + coef(3) = 0.
c                       coef(1) is supposed to be nozero and the maximium 
c                       coefficient is 1 in magnitude
c
c  input arguments:  
c     coef(dp, coef(3),in) - the coeficiects of (ax^2 + bx + c = 0). 
c                            here a = coef(1)
c
c  output arguments:  none
c     numRoots (int, sc, out)  -  number of roots returned. 
c                                     0: no root.
c                                     1: 1 root.
c                                     2: 2 roots.
c     roots (dp, roots(2), out) - the roots.
c
#include "impcom.inc"
c
      double precision coef(3),roots(2),delta
      integer   numRoots

      numRoots = 0

      coef(2) = coef(2)/coef(1)
      coef(3) = coef(3)/coef(1)
      coef(1) = 1.0d0
      delta = coef(2)*coef(2) - 4.0d0*coef(3)
      
      if (abs(delta) .lt. 1.0d-15)then
          numRoots = 1
          roots(1) = -coef(2)*0.5d0
          goto 999
      endif

      if (delta .lt. -1.0d-15) then
          numRoots = 0
          goto 999
      endif

      numRoots = 2
      delta = sqrt(delta)
      roots(1) = 0.5d0*(-coef(2) + delta) 
      roots(2) = 0.5d0*(-coef(2) - delta) 

 999  continue

      return
      end
