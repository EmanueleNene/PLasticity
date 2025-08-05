      subroutine EdgePos (xap,yap,xc,yc,orient)
c *** primary function: this routine determine the relative position of a point (xc,yc) 
c                       with respect of a edge of the polygon.
c
c  input arguments:
c     xc       (dbl,sc,in)          -  x-coordinates of the point 
c     yc       (dbl,sc,in)          -  y-coordinates of the point 
c     xap      (dbl,ar(*),in)       -  x-coordinates of the line end points 
c     yap      (dbl,ar(*),in)       -  y-coordinates of the line end points  
c
c  output arguments:
c     orient   (log,ar(2),out)      -  logical flags, vector of length 2.
c                                      F/F if (xc,yc) is strictly on the right 
c                                      T/F if (xc,yc) is strictly on the left
c                                      F/T if (xc,yc) is on the line segment,
c                                                     and strictly between the end points
c
c
#include "impcom.inc"
c
      logical orient(2)
      double precision  xap(2), yap(2), xc, yc

      double precision  SMALL, SMALLarea, ONEM, ZERO,
     &                  dx1, dy1, dx2, dy2, acr, acl
      parameter (SMALL = 1.0d-10, ONEM = 1.0d0-SMALL, ZERO = 0.0d0)

      dx1 = xap(2) - xap(1)
      dy1 = yap(2) - yap(1)
      dx2 = xc - xap(1)
      dy2 = yc - yap(1)

      SMALLarea = SMALL * (dx1 * dx1 + dy1 * dy1)
      acr = dx1 * dy2 - dx2 * dy1

      if (acr .gt. ZERO) then
          orient(1) = .true.
      else
          orient(1) = .false.
      endif

      orient(2)=.false.

      if (abs(acr).le.SMALLarea) then
          if (abs(dx1).gt.abs(dy1)) then
              acl = dx2 / dx1 
          else
              acl = dy2 / dy1
          endif
          if (acl.gt.SMALL .and. acl.lt.ONEM) orient(2) = .true.
      endif
      
      return
      end

