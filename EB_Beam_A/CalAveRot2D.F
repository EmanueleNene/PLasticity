      subroutine CalAveRot2D (r1,r2,rave)
c --- primary function: compute the average rotation using SLERP
c
c  uses the SLERP procedure (spherical linear interpolation)
c  note, though, we are not using quaternions but rather simple Euler rotations
c
c  ref: http://en.wikipedia.org/wiki/Slerp
c
c  input arguments:
c     r1   (dp,sc,in)     - Euler rotation 1 (theta about n)
c     r2   (dp,sc,in)     - Euler rotation 2
c
c  output arguments:
c     rave (dp,sc,out)    - average rotation
c 
#include "impcom.inc"
c
      double precision  r1,r2,rave,rot1,rot2,q1(0:1),q2(0:1),
     &                  chalf,shalf,angle,qmid(0:1),con,SMALL
      parameter (SMALL = 1.0d-12)
      
c --- convert Euler rotations to quaternions
      rot1 = r1*r1
      if (rot1.lt.SMALL) then         ! very small rotation
          q1(0) = 1.0d0
          q1(1) = 0.0d0
      else
          rot1 = sqrt(rot1)
          chalf = cos(0.5d0*rot1)
          shalf = sin(0.5d0*rot1)
          q1(0) = chalf
          q1(1) = shalf*r1/rot1
      endif

      rot2 = r2*r2
      if (rot2.lt.SMALL) then          ! very small rotation
          if (rot1.lt.SMALL) then      ! catch case when both 
                                       ! rotations are small
              rave = 0.5d0*(r1 + r2)
              goto 999
          endif
          q2(0) = 1.0d0
          q2(1) = 0.0d0
      else
         rot2 = sqrt(rot2)
         chalf = cos(0.5d0*rot2)
         shalf = sin(0.5d0*rot2)
         q2(0) = chalf
         q2(1) = shalf*r2/rot2
      endif

c --- adjust the sign to keep SLERP happy  (q2 and -q2 are the same rotation)
      chalf = q1(0)*q2(0) + q1(1)*q2(1)
      if (chalf.lt.-0.017d0) then                     ! cos(91)
          q2(0) = -q2(0)
          q2(1) = -q2(1)
          chalf = -chalf
      endif

c --- find the angle between the two
      if (chalf.gt.1.0d0) chalf = 1.0d0
      if (chalf.lt.-1.0d0) chalf = -1.0d0
      angle = 2.0d0*acos(chalf)

c --- find the midpoint quaternion using SLERP
      con = 1.0d0/(2.0d0*cos(0.25d0*angle))
      qmid(0) = con*(q1(0) + q2(0))
      qmid(1) = con*(q1(1) + q2(1))

c --- convert quaternion to Euler rotation
      con = qmid(1)*qmid(1)
      if (con.lt.SMALL) then         ! very small rotation
          rave = 0.5d0*(r1 + r2)
      else
          con = sqrt(con)
          chalf = qmid(0)
          if (chalf.gt.1.0d0) chalf = 1.0d0
          if (chalf.lt.-1.0d0) chalf = -1.0d0
          angle = 2.0d0*acos(chalf)
          rave = angle*qmid(1)/con
      endif

 999  return
      end
