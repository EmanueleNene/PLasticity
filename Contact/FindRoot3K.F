      subroutine FindRoot3K (coef, roots, numRoots)
c --- primary function: find all real roots of a degree 3 polynomial
c                       whose leading coefficient is nonzeros and the maximium 
c                       coefficient is 1.0 in magnitude (Kurosh's version)
c
#include "impcom.inc"
#include "constants.inc"
c
      integer   numRoots
      double precision coef(4), roots(3)
      double precision THIRD , QUATER, HALF, TWO, tol1, tol2, tol3
      parameter (QUATER = 0.25d0,
     &           THIRD = 1.0d0/3.0d0, 
     &           HALF  = 0.5d0,      
     &           TWO   = 2.0d0, 
     &           tol1  = 1.0d-10,
     &           tol2  = 2.5d-20,
     &           tol3  = 1.0d-30 )
      
      double precision cm, p,q,r, a,b, del, phi, con1, con2
      integer i

      
      cm = 1.0d0
     &         /max(abs(coef(1)),abs(coef(2)),abs(coef(3)),abs(coef(4)))
      
      numRoots = 0
      
      if (abs(coef(1) *cm ) .le. tol1) then
          if (abs(coef(2) *cm) .le. tol1) then
              numRoots = 1
              roots(1) = -coef(4)/coef(3)
              goto 999
          else
              coef(2) = coef(2) * cm
              coef(3) = coef(3) * cm 
              coef(4) = coef(4) * cm
              call FindRoot2 (coef(2),roots(1),numRoots)
              goto 999
          endif
      endif
      
      p = coef(2)/coef(1)
      q = coef(3)/coef(1)
      r = coef(4)/coef(1)
      a = q - THIRD*p*p
      b = THIRD*THIRD*THIRD * (TWO*p*p*p - 9.0d0*p*q + 27.0d0*r)
      
      del = QUATER*b*b + a*a*a*THIRD*THIRD*THIRD
      
      if (del .gt. tol1) then
          numroots = 1
          con1 = -HALF*b + sqrt(del)
          con2 = -HALF*b - sqrt(del)
          roots(1) = sign(exp(THIRD*log(abs(con1))),con1) +
     &               sign(exp(THIRD*log(abs(con2))),con2) - 
     &               THIRD*p
      elseif (del .gt. -tol2) then
          if (b .gt. tol3) then
              numroots = 2
              roots(1) = -TWO*sqrt(-THIRD*a) - THIRD*p
              roots(2) = sqrt(-THIRD*a) - THIRD*p
          elseif (b .gt. -tol3) then
              numroots = 1
              roots(1) = -THIRD*p
          else
              numroots = 2
              roots(1) = TWO * sqrt(-THIRD*a) - THIRD*p
              roots(2) = -sqrt(-THIRD*a) - THIRD*p
          endif
      else
          numroots = 3
          if (del .gt. 0.0d0) then
            if (b .gt. 0.0d0) then
              phi = acos(-sqrt((QUATER*b*b)/(a*a*a*THIRD*THIRD*THIRD)))
            else
              phi = acos(sqrt((QUATER*b*b)/(a*a*a*THIRD*THIRD*THIRD)))
            endif
          else
            if (b .gt. 0.0d0) then
              phi = acos(-sqrt((QUATER*b*b)/(-a*a*a*THIRD*THIRD*THIRD)))
            else
              phi = acos(sqrt((QUATER*b*b)/(-a*a*a*THIRD*THIRD*THIRD)))
            endif
          endif
          do i = 1,3
             roots(i) = TWO*sqrt(-a*THIRD)
     &                  *cos(THIRD*(phi + TWO*dble(i-1)*W_Pi)) - THIRD*p
          enddo
          
c          if ( abs(roots(2)-roots(3)) .le. tol0 ) then
c              numRoots = 2
c              roots(2) = HALF * ( roots(2) + roots(3) )
c          endif
          
      endif
C
 999  continue
 
      return
      end