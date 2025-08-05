      function CalDist2P (a,b)
c --- function: find the distance squared between two points
c
c  input arguments:
c     a        (dp,ar(3),in)     - x-y-z position of first  point
c     b        (dp,ar(3),in)     - x-y-z position of second point
c
c  output arguments:
c     CalDist2P   (dp,sc,out)       - distance squared between the two points
c
#include "impcom.inc"
c
      double precision   a(3),b(3),CalDist2P
c

      CalDist2P = sqrt((a(1)-b(1))**2 + (a(2)-b(2))**2 + (a(3)-b(3))**2)

      return
      end
