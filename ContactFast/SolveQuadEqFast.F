      subroutine SolveQuadEqFast (a,b,u,v,f1,pf1_s,pf1_t,f2,pf2_s,pf2_t)
c --- primary function: solve system of quadratic equations
c                       coef(1)*x^2 + coef(2)*x + coef(3) = 0.
c
c  input arguments:
c     a  (dp,ar(4,3),in) - the coeficiects for the first beam
c     b  (dp,ar(4,3),in) - the coeficiects for the second beam
c
#include "impcom.inc"
c
      double precision a(4,3),b(4,3),u,v,f1,pf1_s,pf1_t,f2,pf2_s,pf2_t
      double precision d1,d2,d3,u2,v2,u3,v3,u1v1,u2v1,u1v2,
     &                 c3,c4,c5,c10,c11,c13,c14,c15
c

      u2 = u*u
      v2 = v*v
      u3 = u*u2     
      v3 = v*v2
      u1v1 = u*v
      u2v1 = u2*v
      u1v2 = u*v2

      d1 = a(1,1)-b(1,1)
      d2 = a(1,2)-b(1,2)
      d3 = a(1,3)-b(1,3)

      c3 = 2.d0*(a(3,1)*a(3,1)+a(3,2)*a(3,2)+a(3,3)*a(3,3))
      c4 = 3.d0*(a(2,1)*a(3,1)+a(2,2)*a(3,2)+a(2,3)*a(3,3))
      c5 = a(2,1)*a(2,1)+a(2,2)*a(2,2)+a(2,3)*a(2,3)+
     &     2.d0*(a(3,1)*d1+a(3,2)*d2+a(3,3)*d3)

      c10= -2.d0*(a(3,1)*b(3,1)+a(3,2)*b(3,2)+a(3,3)*b(3,3))
      c11= -(a(2,1)*b(3,1)+a(2,2)*b(3,2)+a(2,3)*b(3,3))
      c13= -2.d0*(a(3,1)*b(2,1)+a(3,2)*b(2,2)+a(3,3)*b(2,3))
      c14= -(a(2,1)*b(2,1)+a(2,2)*b(2,2)+a(2,3)*b(2,3))
      c15 = a(2,1)*d1+a(2,2)*d2+a(2,3)*d3

c --- calculate residuals and slopes
      f1 = c3*u3+c4*u2+c5*u+c10*u1v2+c11*v2+c13*u1v1+c14*v+c15
      pf1_s = 3.d0*c3*u2+2.d0*c4*u+c5+c10*v2+c13*v
      pf1_t = 2.d0*c10*u1v1+2.d0*c11*v+c13*u+c14

      c3 = 2.d0*(b(3,1)*b(3,1)+b(3,2)*b(3,2)+b(3,3)*b(3,3))
      c4 = 3.d0*(b(2,1)*b(3,1)+b(2,2)*b(3,2)+b(2,3)*b(3,3))
      c5 = b(2,1)*b(2,1)+b(2,2)*b(2,2)+b(2,3)*b(2,3) -
     &     2.d0*(b(3,1)*d1+b(3,2)*d2+b(3,3)*d3)

      c11= -(b(2,1)*a(3,1)+b(2,2)*a(3,2)+b(2,3)*a(3,3))
      c13= -2.d0*(b(3,1)*a(2,1)+b(3,2)*a(2,2)+b(3,3)*a(2,3))
      c15 = -(b(2,1)*d1+b(2,2)*d2+b(2,3)*d3)

c --- calculate residuals and slopes
      f2 = c3*v3+c4*v2+c5*v+c10*u2v1+c11*u2+c13*u1v1+c14*u+c15
      pf2_s = 2.d0*c10*u1v1+2.d0*c11*u+c13*v+c14
      pf2_t = 3.d0*c3*v2+2.d0*c4*v+c5+c10*u2+c13*u

      return
      end

