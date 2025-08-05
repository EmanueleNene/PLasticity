      subroutine BndPts(W1,W2,Ang,XY)
c --- primary function: Finds the corner points of the bonded area
c
c input arguments:
c   W1      (dp,sc,in)  -Width of beam 1
c   W2      (dp,sc,in)  -Width of beam 2
c   Ang     (dp,sc,in)  -Angle between two beams
c
c output arguments:
c   XY      (dp,arr(4,2),out) -Coordinates of four corner points
c
#include "impcom.inc"
c
      double precision    W1,W2,Ang
      double precision    XY(4,2)
      
      double precision    a,c,w
      double precision    x1, x2
      
      a = 0.5d0/sin(ang)
      c = cos(Ang)
      w = 0.5d0*W1
      
      x1 = a*(W2 - W1*c)
      x2 = a*(W2 + W1*c)
      
          XY(1,1) = -x2
          XY(2,1) =  x1
          XY(3,1) =  x2
          XY(4,1) = -x1
      
          XY(1,2) = -w
          XY(2,2) = -w
          XY(3,2) =  w
          XY(4,2) =  w
      
      return
      end