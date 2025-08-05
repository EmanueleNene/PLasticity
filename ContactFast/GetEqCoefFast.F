      subroutine GetEqCoefFast (xa,ya,za,xb,yb,zb,a,b)
c --- primary function: get coeffs for system of quadratic equations
c                       coef(1)*x^2 + coef(2)*x + coef(3) = 0.
c
#include "impcom.inc"
c
      double precision xa(3),ya(3),za(3),xb(3),yb(3),zb(3),
     &                 La,Lb
      double precision a(4,3),b(4,3)
      double precision ZERO,HALF
      parameter (ZERO=0.d0,HALF=0.5d0)

      La = (xa(1)-xa(2))**2+(ya(1)-ya(2))**2+(za(1)-za(2))**2
      Lb = (xb(2)-xb(3))**2+(yb(2)-yb(3))**2+(zb(2)-zb(3))**2
      if (La.lt.Lb) La = Lb
      La = sqrt(La)


      a(1,1) = xa(3)/La
      a(2,1) = HALF*(xa(2)-xa(1))/La
      a(3,1) = HALF*(xa(1)+xa(2)-2.d0*xa(3))/La
      a(4,1) = ZERO
            
      a(1,2) = ya(3)/La
      a(2,2) = HALF*(ya(2)-ya(1))/La
      a(3,2) = HALF*(ya(1)+ya(2)-2.d0*ya(3))/La      
      a(4,2) = ZERO
      
      a(1,3) = za(3)/La
      a(2,3) = HALF*(za(2)-za(1))/La
      a(3,3) = HALF*(za(1)+za(2)-2.d0*za(3))/La
      a(4,3) = ZERO


      b(1,1) = xb(3)/La
      b(2,1) = HALF*(xb(2)-xb(1))/La
      b(3,1) = HALF*(xb(1)+xb(2)-2.d0*xb(3))/La
      b(4,1) = ZERO

      b(1,2) = yb(3)/La
      b(2,2) = HALF*(yb(2)-yb(1))/La
      b(3,2) = HALF*(yb(1)+yb(2)-2.d0*yb(3))/La
      b(4,2) = ZERO

      b(1,3) = zb(3)/La
      b(2,3) = HALF*(zb(2)-zb(1))/La
      b(3,3) = HALF*(zb(1)+zb(2)-2.d0*zb(3))/La
      b(4,3) = ZERO

      return
      end
