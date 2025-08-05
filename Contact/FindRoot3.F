      subroutine FindRoot3 (coef, roots, numRoots)
c --- primary function: find all real roots of a degree 3 polynomial
c                       whose leading coefficient is nonzeros and the maximium 
c                       coefficient is 1.0 in magnitude
c
#include "impcom.inc"
c
      double precision coef(4),s,t,b,c,aa1,d,aa2,roots(3),xmax
      integer   numRoots

      numRoots = 0
      xmax = max(abs(coef(1)),abs(coef(2)),abs(coef(3)),abs(coef(4)))
      coef(1) = coef(1)/xmax
      coef(2) = coef(2)/xmax   
      coef(3) = coef(3)/xmax   
      coef(4) = coef(4)/xmax   

      if (abs(coef(1)).le.1.0d-10.and.abs(coef(2)).le.1.0d-10) then
          numRoots = 1
          roots(1) = -coef(4)/coef(3)
          goto 999
      endif

      if (abs(coef(1)).le.1.0d-10) then
          call FindRoot2 (coef(2),roots(1),numRoots)
          goto 999
      endif

      coef(2) = coef(2)/coef(1)
      coef(3) = coef(3)/coef(1)
      coef(4) = coef(4)/coef(1)
      coef(1) = 1.0d0
   
      s = coef(2)/3.0d0
      t = s*coef(2)
      b = 0.5d0*(s*(t/1.5d0 - coef(3)) + coef(4))
      t = (t-coef(3))/3.0d0
      if (abs(t).gt.1.0d+14) then
          roots(1) = 0
          numRoots = 0
          goto 999
      endif 
      c = t*t*t
      d = b*b - c
      if (d.lt.0.0d0) goto 107
      aa1 = sqrt(d) + abs(b)
      d = 0.0d0
      if (aa1.gt.1.d-8) then
          aa2 = 1.0d0/3.d0*log(aa1)
          d = exp(aa2)
      endif
      if (d.eq.0.0) goto 105
      if (b.le.0.0d0) goto 103
      b = -d
      goto 104
 103  b = d
 104  c = t/b
 105  continue
      d = sqrt(0.75d0)*(b-c)
      b = b + c
      c = -0.5d0*b-s
      roots(1) = b - s
      numRoots = 1      
      if (abs(d).le.1.d-8*abs(c)) then
          numRoots = 2      
          roots(2) = c
      endif
      goto 999
 107  continue
c
      if(b.ne.0.d0) goto 108
      d = atan(1.d0)/1.5d0
      goto 109
 108  d = atan(sqrt(-d)/abs(b))/3.d0
 109  if (b.ge.0.0d0) goto 110
      b = sqrt(t)*2.d0
      goto 111
 110  b = sqrt(t)*(-2.d0)
 111  c = cos(d)*b
      t = -sqrt(0.75d0)*sin(d)*b - 0.5d0*c
      d = -t-c-s
      c = c-s
      t = t-s
      if (abs(c).le.abs(t)) goto 112
      roots(3) = c
      goto 113
 112  roots(3) = t
      t = c
 113  if(abs(d).le.abs(t)) goto 114
      roots(2) = d
      goto 115
 114  roots(2) = t
      t = d
 115  roots(1) = t
      numRoots = 3      
C
 999  continue
 
      return
      end
