      function IfContain (nPoly,xyPoly,xy,option)
c *** primary function: determine if a point (xc,yc) is in inside a polygon
c                       by Darel Rex Finley's method
c
c  input arguments:
c     nPoly       (int,sc,in)        - number of segments
c     xyPoly      (dp,ar(2,*),in)    - xy values of start point
c     xy          (dp,ar(2),in)      - xy values of the point
c     option      (int,sc,in)        - option of treating the edge node (0 - inside, 1 - outside)
c                                            
c  output arguments:
c     IfContain  (log,sc,out)        - .true.  = inside
c                                      .false. = outside
c                                             
#include "impcom.inc"
c
      logical   IfContain
      integer   nPoly, option
      double precision
     &          xyPoly(2,nPoly),xy(2)
 
       integer   i
      double precision
     &          x1, y1, x2, y2, dx1, dy1, dx2, dy2, 
     &          acr, acl, SMALLarea
      logical   inside

      double precision SMALL, ONEM
      parameter (SMALL = 1.0d-10, ONEM = 1.0d0-SMALL)

      inside = .false.
      
      x1 = xyPoly(1,nPoly)
      y1 = xyPoly(2,nPoly)

      do i = 1,nPoly

c --- check if the point coincide with one vertex of the polygon
         if (((x1 - xy(1))**2+(y1 - xy(2))**2).lt.SMALL) then
             if (option.eq.1) then
                 inside = .false.
             else
                 inside = .true.
             endif    
             goto 999
         endif

         x2 = xyPoly(1,i)
         y2 = xyPoly(2,i)

c --- check if the point is on the edge
         dx1 = x2 - x1
         dy1 = y2 - y1
         dx2 = xy(1) - x1
         dy2 = xy(2) - y1

         SMALLarea = SMALL * (dx1 * dx1 + dy1 * dy1)
         acr = dx1 * dy2 - dx2 * dy1

         if (abs(acr).le.SMALLarea) then
             if (abs(dx1).gt.abs(dy1)) then
                 acl = dx2 / dx1 
             else
                 acl = dy2 / dy1
             endif
             if (acl.gt.SMALL .and. acl.lt.ONEM) then
                 if (option.eq.1) then
                     inside = .false.
                 else
                     inside = .true.
                 endif    
                 goto 999
             endif
         endif

         if (y2.lt.xy(2) .and. y1.ge.xy(2) .or. 
     &       y1.lt.xy(2) .and. y2.ge.xy(2) .and. 
     &      (x2.le.xy(1) .or.  x1.le.xy(1))) then
             if (x2+(xy(2)-y2)/(y1-y2)*(x1-x2).lt.xy(1)) then
                 inside = .not.inside
             endif
         endif
         x1 = x2
         y1 = y2
      enddo


 999  continue
       
      IfContain = inside

      return
      end      
