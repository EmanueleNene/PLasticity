      subroutine NodeAtDist(len,C,xyz)
c --- primary function: Finds the xyz coordinates of a point at length=len on a curve defined with C coefficients
c     using newton-raphson iterations.
c     Input:
c       C           Segment Coefficients (x = C(1) + C(2)*t + C(3)*t^2 + C(4)*t^3, y = ... , z = ...)
c       len         length of travel from beginning node to desired point
c     OutPut:
c       xyz         xyz coordinates of required point
c
#include "impcom.inc"
#include "constants.inc"
c
      double precision TOL
      parameter (TOL = 1.0d-7)
      
      double precision C(12), xyz(3), len
      double precision SegLength

      double precision l,lp, xyzp(3), t, t0, t1
      integer   count

      t0 = 0.0d0
      t1 = 1.0d0
      t  = len / SegLength(C,t0,t1)
      count = 0
      t1 = t + TOL + 1.0d0
      
      do while ( (abs(t1-t) .gt. TOL) .and. (count .lt. 20) )
         count = count + 1
         xyzp(1) = C(2 ) + 2.0d0*C(3 )*t + 3.0d0*C(4 )*t*t
         xyzp(2) = C(6 ) + 2.0d0*C(7 )*t + 3.0d0*C(8 )*t*t
         xyzp(3) = C(10) + 2.0d0*C(11)*t + 3.0d0*C(12)*t*t
         l  = SegLength(C,t0,t)
         lp = (xyzp(1)*xyzp(1)+ xyzp(2)*xyzp(2)+ xyzp(3)*xyzp(3))**0.5d0
         t1 = t
         t = t - (l-len)/lp
      enddo
      if ( (t .lt. 0) .or. (t .gt. 1) ) t = len / SegLength(C,t0,t1)

      xyz(1) = C(1) +  C(2)*t  + C(3)*t*t  + C(4)*t*t*t
      xyz(2) = C(5) +  C(6)*t  + C(7)*t*t  + C(8)*t*t*t
      xyz(3) = C(9) + C(10)*t + C(11)*t*t + C(12)*t*t*t
      
      return
      end