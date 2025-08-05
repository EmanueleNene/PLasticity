      subroutine GetRotMat (phi, r)
c -- primary function:  get rotation matrix from rotation angles
c
c  input argument:
c        phi      (dp,ar(3),in)     - rotation angles
c
c  output argument:
c        r        (dp,ar(3,3),out)  - rotation matrix
c
#include "impcom.inc"
c
      double precision phi(3), r(3,3), p(3), theta, cost, sint

      theta = phi(1)*phi(1) + phi(2)*phi(2) + phi(3)*phi(3)

      if (theta .gt. 1.0d-24) then
          theta = sqrt(theta)
          p(1) = phi(1)/theta
          p(2) = phi(2)/theta
          p(3) = phi(3)/theta
          cost = cos(theta)
          sint = sin(theta)
          r(1,1) = cost + (1.0d0-cost)*p(1)*p(1)
          r(1,2) = (1.0d0-cost)*p(1)*p(2) - sint*p(3)
          r(1,3) = (1.0d0-cost)*p(1)*p(3) + sint*p(2)
          r(2,1) = (1.0d0-cost)*p(2)*p(1) + sint*p(3)
          r(2,2) = cost + (1.0d0-cost)*p(2)*p(2) 
          r(2,3) = (1.0d0-cost)*p(2)*p(3) - sint*p(1)
          r(3,1) = (1.0d0-cost)*p(3)*p(1) - sint*p(2)
          r(3,2) = (1.0d0-cost)*p(3)*p(2) + sint*p(1)
          r(3,3) = cost + (1.0d0-cost)*p(3)*p(3)
      else
c --- no rotation
          call v_zero(r(1,1),9)
          r(1,1) = 1.0d0
          r(2,2) = 1.0d0
          r(3,3) = 1.0d0
      endif

      return
      end
