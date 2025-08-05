      subroutine CalSect_6(intp, nnodi, parInt)
c --- primary function: calculate section data
c
#include "impcom.inc"
c
      integer MNODI
      parameter(MNODI = 9)
      integer intp, nnodi
      double precision parInt(MNODI,3),g,h,ums,umt,ups,upt,umss,umtt
      double precision ONE, HALF, TWO, gh(2,4),
     &     g1, g2, g3, h1, h2, h3, d1, d2, d3, FOURTH
      parameter (ONE=1.0d0, HALF=0.5D0, TWO=2.0d0,
     &     FOURTH=0.25D0)
      data gh/
     &     -0.577350269189626d0,-0.577350269189626d0,
     &      0.577350269189626d0,-0.577350269189626d0,
     &      0.577350269189626d0, 0.577350269189626d0,
     &     -0.577350269189626d0, 0.577350269189626d0/
c

      g = gh(1,intp)
      h = gh(2,intp)
      if (nnodi.eq.9) then
          g1 = HALF*g*(g-ONE)
          g2 = HALF*g*(g+ONE)
          g3 = ONE-g*g
          h1 = HALF*h*(h-ONE)
          h2 = HALF*h*(h+ONE)
          h3 = ONE-h*h
          parInt(1,1) = g1*h1
          parInt(2,1) = g2*h1
          parInt(3,1) = g2*h2
          parInt(4,1) = g1*h2
          parInt(5,1) = g3*h1
          parInt(6,1) = g2*h3
          parInt(7,1) = g3*h2
          parInt(8,1) = g1*h3
          parInt(9,1) = g3*h3
          d1 = g-HALF
          d2 = g+HALF
          d3 = -TWO*g
          parInt(1,2) = d1*h1
          parInt(2,2) = d2*h1
          parInt(3,2) = d2*h2
          parInt(4,2) = d1*h2
          parInt(5,2) = d3*h1
          parInt(6,2) = d2*h3
          parInt(7,2) = d3*h2
          parInt(8,2) = d1*h3
          parInt(9,2) = d3*h3
          d1 = h-HALF
          d2 = h+HALF
          d3 = -TWO*h
          parInt(1,3) = g1*d1
          parInt(2,3) = g2*d1
          parInt(3,3) = g2*d2
          parInt(4,3) = g1*d2
          parInt(5,3) = g3*d1
          parInt(6,3) = g2*d3
          parInt(7,3) = g3*d2
          parInt(8,3) = g1*d3
          parInt(9,3) = g3*d3
      elseif (nnodi.eq. 8) then
          ums = ONE - g
          umt = ONE - h
          ups = ONE + g
          upt = ONE + h
          umss = ums*ups
          umtt = umt*upt
          parInt(1,1) = FOURTH*ums*umt*(-g-h-ONE)
          parInt(2,1) = FOURTH*ups*umt*( g-h-ONE)
          parInt(3,1) = FOURTH*ups*upt*( g+h-ONE)
          parInt(4,1) = FOURTH*ums*upt*(-g+h-ONE)
          parInt(5,1) = HALF*umss*umt
          parInt(6,1) = HALF*ups*umtt
          parInt(7,1) = HALF*umss*upt
          parInt(8,1) = HALF*ums*umtt
          parInt(1,2) = FOURTH*umt*(TWO*g + h)
          parInt(2,2) = FOURTH*umt*(TWO*g - h)
          parInt(3,2) = FOURTH*upt*(TWO*g + h)
          parInt(4,2) = FOURTH*upt*(TWO*g - h)
          parInt(5,2) = -g*umt
          parInt(6,2) = HALF*umtt
          parInt(7,2) = -g*upt
          parInt(8,2) = -HALF*umtt
          parInt(1,3) = FOURTH*ums*(TWO*h + g)
          parInt(2,3) = FOURTH*ups*(TWO*h - g)
          parInt(3,3) = FOURTH*ups*(TWO*h + g)
          parInt(4,3) = FOURTH*ums*(TWO*h - g)
          parInt(5,3) = -HALF*umss
          parInt(6,3) = -h*ups
          parInt(7,3) = HALF*umss
          parInt(8,3) = -h*ums
      endif

      return
      end
