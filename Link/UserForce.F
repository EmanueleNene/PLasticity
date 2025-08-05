      subroutine UserForce (Disp,c,Stiff,Force,key)
c --- primary function: define displacement-force relations and the tangent stiffness
c     
c  input arguments:
c     Disp  (dp,sc,in)        - calculated displacement
c     c     (dp,ar(5),in)     - constants
c     key   (int,sc,in)       - keyopt(3) on the spring element
c     
c  output arguments:
c     Stiff  (dp,sc,inout)    - element number
c     Force  (dp,sc,inout)    - array of element data
c
#include "impcom.inc"
#include "constants.inc"
c     
      integer key
      double precision  Disp, Stiff, Force, c(3)
      
      if (key.eq.1) then
          Force = c(1)*sinh(c(2)*Disp)
          Stiff = c(2)*c(1)*cosh(c(2)*Disp)
      endif
      
      return
      end


