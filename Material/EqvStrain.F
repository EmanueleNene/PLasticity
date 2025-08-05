      function EqvStrain (ep,nu)
c --- primary function: combines strain components to give an "overall" strain
c                       used in creep and plasticity calculations
c
c  input arguments:
c     ep       (dp,ar(6),in)       - the strain components
c     nu       (dp,sc,in)          - Poisson's ratio
c
c  output arguments:
c     EqvStrain     (dp,func,out)  - the combined strain value
c
c     Expression in general case is:
c
c         ___________________________________________________________________
c        /1             2             2             2  1    2      2      2
c    \  / -*((ep  - ep ) + (ep  - ep ) + (ep  - ep ) + -*(ep   + ep   + ep ))
c     \/  2     1     2       2     3       3     1    2   4      5      6
c    ------------------------------------------------------------------------
c                                (1 + nu)
c
c    NOTE: no tangent stresses
c
#include "impcom.inc"
c
      double precision  ep(6),nu,EqvStrain

      EqvStrain = 2*ep(1)**2

      if (EqvStrain.gt.0.0d0)
     &    EqvStrain = sqrt (0.5d0*EqvStrain)/(1.0d0+nu)
c
      return
      end
