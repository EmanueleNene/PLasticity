      function CheckJacBeam (xyz)
c --- primary function: compute jacobian ratio for 3-node beam elements
c
c  input arguments:
c     xyz    (dp,ar(3,3),in)   - node coordinates (3 coordinates (rows) x 3 nodes)
c
c  output arguments:
c     ratj   (dp,sc,out)       - ratio of jacobians evaluated at nodes
c                                if -1.d0 => NOT ACCEPTABLE (zero or negative)
c
c  NOTE: not pasing the test does not mean the element produces totaly wrong results
c        only the accuracy is affected
c
#include "impcom.inc"
c
      integer  i,k
      double precision  CheckJacBeam,ratj,jvect(3,3),xyz(3,3),
     &                  sdln3nodes(3,3),con,dot12,dot23,dot31

c --- shape functions for 3-node beam (3 integration points per length)
c     IP 1:  S = -1.0000000000000
c     IP 2:  S =  0.0000000000000
c     IP 3:  S =  1.0000000000000
      data  ((sdln3nodes(i,k),i=1,3),k=1,3) /
     & -1.5000000000000D0,  -0.5000000000000D0,   2.0000000000000D0,
     & -0.5000000000000D0,   0.5000000000000D0,   0.0000000000000D0,
     &  0.5000000000000D0,   1.5000000000000D0,  -2.0000000000000D0/


#if defined(DEBUG)
      call SubBeg('CheckJacBeam')
#endif

      i = 0

      jvect(1,1) = sdln3nodes(1,1)*xyz(1,1)
     &           + sdln3nodes(2,1)*xyz(1,2)
     &           + sdln3nodes(3,1)*xyz(1,3)
      jvect(2,1) = sdln3nodes(1,1)*xyz(2,1)
     &           + sdln3nodes(2,1)*xyz(2,2)
     &           + sdln3nodes(3,1)*xyz(2,3)
      jvect(3,1) = sdln3nodes(1,1)*xyz(3,1)
     &           + sdln3nodes(2,1)*xyz(3,2)
     &           + sdln3nodes(3,1)*xyz(3,3)

      jvect(1,2) = sdln3nodes(1,2)*xyz(1,1)
     &           + sdln3nodes(2,2)*xyz(1,2)
     &           + sdln3nodes(3,2)*xyz(1,3)
      jvect(2,2) = sdln3nodes(1,2)*xyz(2,1)
     &           + sdln3nodes(2,2)*xyz(2,2)
     &           + sdln3nodes(3,2)*xyz(2,3)
      jvect(3,2) = sdln3nodes(1,2)*xyz(3,1)
     &           + sdln3nodes(2,2)*xyz(3,2)
     &           + sdln3nodes(3,2)*xyz(3,3)

      jvect(1,3) = sdln3nodes(1,3)*xyz(1,1)
     &           + sdln3nodes(2,3)*xyz(1,2)
     &           + sdln3nodes(3,3)*xyz(1,3)
      jvect(2,3) = sdln3nodes(1,3)*xyz(2,1)
     &           + sdln3nodes(2,3)*xyz(2,2)
     &           + sdln3nodes(3,3)*xyz(2,3)
      jvect(3,3) = sdln3nodes(1,3)*xyz(3,1)
     &           + sdln3nodes(2,3)*xyz(3,2)
     &           + sdln3nodes(3,3)*xyz(3,3)

      dot12 = jvect(1,1)*jvect(1,2)
     &      + jvect(2,1)*jvect(2,2)
     &      + jvect(3,1)*jvect(3,2)
      dot23 = jvect(1,2)*jvect(1,3)
     &      + jvect(2,2)*jvect(2,3)
     &      + jvect(3,2)*jvect(3,3)
      dot31 = jvect(1,3)*jvect(1,1)
     &      + jvect(2,3)*jvect(2,1)
     &      + jvect(3,3)*jvect(3,1)

      if (dot12.gt.0.0d0 .and. dot23.gt.0.0d0 .and. dot31.gt.0.0d0) then
          ratj = dot12/dot23
          if (ratj .lt. 1.0d0)  ratj = 1.0d0/ratj
          con = dot23/dot31
          if (con .lt. 1.0d0)  con = 1.0d0/con
          if (con .gt. ratj)   ratj = con
          con = dot31/dot12
          if (con .lt. 1.0d0)  con = 1.0d0/con
          if (con .gt. ratj)   ratj = con
      else
          ratj = -1.0d0
      endif
      CheckJacBeam = ratj

#if defined(DEBUG)
      call SubEnd('CheckJacBeam')
#endif

      return
      end
