      function SegLength(C,t0,t1)
c --- primary function: calculate length of (part of a / a whole) segment using Gauss quadrature method
c
#include "impcom.inc"
#include "constants.inc"
c
      double precision  C(12), t0, t1, SegLength
      
      integer i
      double precision  t(5),w(5),a,b

      a = (t1-t0) * 0.5d0
      b = (t1+t0) * 0.5d0
      
      t(1) = 0.00000000000000000000000000000000d0 * a + b
      t(2) = 0.53846931010568309103631442070021d0 * a + b
      t(3) =-0.53846931010568309103631442070021d0 * a + b
      t(4) = 0.90617984593866399279762687829939d0 * a + b
      t(5) =-0.90617984593866399279762687829939d0 * a + b
      w(1) = 0.56888888888888888888888888888889d0
      w(2) = 0.47862867049936646804129151483564d0
      w(3) = 0.47862867049936646804129151483564d0
      w(4) = 0.23692688505618908751426404071992d0
      w(5) = 0.23692688505618908751426404071992d0
      
      SegLength = 0.0d0
      do i=1,5
         SegLength = SegLength +  w(i) * (
     &           (C(2) +  2.0d0*C(3) *t(i) + 3.0d0*C(4) *t(i)*t(i))**2 +
     &           (C(6) +  2.0d0*C(7) *t(i) + 3.0d0*C(8) *t(i)*t(i))**2 +
     &           (C(10) + 2.0d0*C(11)*t(i) + 3.0d0*C(12)*t(i)*t(i))**2
     &                                                          )**0.5d0
      enddo
      SegLength = SegLength * a

      return
      end