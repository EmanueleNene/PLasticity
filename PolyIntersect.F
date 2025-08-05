      function PolyIntersect(Cf,xyPoly,nPoly,xyz)
c --- primary function: Finds the coordinates of the intersection point between the segmentand the polygon
c     Input:
c       Cf              Segment Coefficients
c       xyPoly          coordinates of polygon keypoints
c       nPoly           Number of polygon sides
c     OutPut:
c       xyz             coordinates of intersection point
c       PolyIntersect   neutral parametrization value at intersection point
c
#include "impcom.inc"
#include "constants.inc"
c
      double precision   PolyIntersect, TOL
      integer   maxIter
      parameter( TOL        = 1.0d-8,
     &           maxIter    = 20 )

      integer nPoly
      double precision Cf(12), xyPoly(2,nPoly), xyz(3)
      
      double precision  x1,y1,z1,x2,y2,z2
      double precision  m, mp, x
      double precision  a, b, c, d, f, f1, t, t0, tt
      integer   i, j, count
      
      x1 = Cf(1)
      x2 = Cf(1) + Cf(2) + Cf(3) + Cf(4)
      y1 = Cf(5)
      y2 = Cf(5) + Cf(6) + Cf(7) + Cf(8)
      z1 = Cf(9)
      z2 = Cf(9) + Cf(10)+ Cf(11)+ Cf(12)
      
      m = (y2-y1)/(x2-x1)
      
      tt = 2.0d0
      t = -1.0d0
      do i=1,nPoly
         j = mod(i,nPoly) +1
         if (xyPoly(1,j) .eq. xyPoly(1,i)) then
             x = xyPoly(1,i)
         else
             mp = (xyPoly(2,j)-xyPoly(2,i))/(xyPoly(1,j)-xyPoly(1,i))
             x  = (m*x1 - mp*xyPoly(1,i)+ xyPoly(2,i) - y1)/(m-mp)
         endif
         if ( ( x .ge. min(xyPoly(1,i),xyPoly(1,j)) ) .and.
     &        ( x .le. max(xyPoly(1,i),xyPoly(1,j)) ) .and.
     &        ( x .ge. min(x1,x2) ) .and.
     &        ( x .le. max(x1,x2) ) ) then
           if (xyPoly(1,j) .eq. xyPoly(1,i)) then
               a = Cf(1) - xyPoly(1,i)
               b = Cf(2)
               c = Cf(3)
               d = Cf(4)
           else
               a = Cf(5) - mp*Cf(1) + mp*xyPoly(1,i) - xyPoly(2,i)
               b = Cf(6) - mp*Cf(2)
               c = Cf(7) - mp*Cf(3)
               d = Cf(8) - mp*Cf(4)
           endif
          
           count = 1
           t = 0.5d0
           t0 = t + TOL + 1.0d0
           do while ( (abs(t-t0) .gt. TOL) .and. (count .lt. maxIter) )
              t0 = t
              f = a + b*t + c*t*t + d*t*t*t
              f1 = b   + 2.0d0*c*t + 3.0d0*d*t*t
              t = t - f/f1
              count = count +1
           enddo
           if ( (t .gt. 0) .and. (t .lt. 1) ) tt = min(tt,t)
        endif
      enddo
      t = tt
      if ( (t .le. 0) .or. (t .ge. 1) ) t=0.5d0
100   continue
      
      xyz(1) = Cf(1) +  Cf(2)*t  + Cf(3)*t*t  + Cf(4)*t*t*t
      xyz(2) = Cf(5) +  Cf(6)*t  + Cf(7)*t*t  + Cf(8)*t*t*t
      xyz(3) = Cf(9) + Cf(10)*t + Cf(11)*t*t + Cf(12)*t*t*t

      PolyIntersect = t
      return
      end