      function IfInside (nPoly,xyPoly,xc,yc)
c *** primary function: determine if a point (xc,yc) is in inside a polygon
c
c  input arguments:
c     nPoly    (int,sc,in)          -  number of nodes of the polygon 
c     xyPoly   (db,ar(2,*),in)      -  x,y-coordinates of the polygon
c     xc       (db,ar(*),in)        -  x-coordinates of the point 
c     yc       (db,ar(*),in)        -  y-coordinates of the point 
c
c  output arguments:
c     IfInside (log,sc,fun)         -  .true. - inside; .false. - outside
c
c  local arguments:
c     locp     (log,ar(3),out)      -  logical flags, vector of length 3:
c                                      F/F/F if (xc,yc) is strictly outside
c                                      T/F/F if (xc,yc) is strictly inside
c                                      F/T/F if (xc,yc) is on the edge
c                                      F/F/T if (xc,yc) coincides with a vertex
c
c
#include "impcom.inc"
c
      logical IfInside, locp(3)
      integer nPoly
      double precision  xyPoly(2,nPoly), xc, yc

      logical res(2)
      integer i, j, in
      double precision SMALL, xaP(2),yaP(2), dx, dy, ds
      parameter (SMALL = 1.0d-20)

      locp(1) = .true.
      IfInside = .false.
c --- loop over edges of the polygon

      do 200 i = 1,nPoly
         in = mod(i,nPoly) + 1
         locp(2) = .false.
         locp(3) = .false.

         xaP(1) = xyPoly(1,i)
         xaP(2) = xyPoly(1,in)
         yaP(1) = xyPoly(2,i)
         yaP(2) = xyPoly(2,in)

         call EdgePos (xaP,yaP,xc,yc,res)
         locp(1) = locp(1).and.res(1)
         
c --- if not in the counter clock wise direction
         if (.not.(locp(1))) then
             locp(2) = res(2)

c --- if the point is on the edge
             if (locp(2)) goto 999

c --- check if the point coincide with one vertex of the polygon
             do 205 j = 1,2
                dx = xaP(j) - xc
                dy = yaP(j) - yc
                ds = dx * dx + dy * dy
                if (ds.lt.SMALL) then
                    locp(3) = .true.
                    goto 999
                endif
 205         continue
         endif
 200  continue

 999  continue
      
      IfInside = locp(1).or.locp(2).or.locp(3)
           
      return
      end

