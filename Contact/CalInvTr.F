      subroutine CalInvTr (nfacn,x,y,z,ks,et,vn,tr)
c --- primary function: determine invariant of transformation matrix 
c
c  input arguments:
c     nfacn    (int,sc,in)       - number of nodes
c     xyz      (dp,ar(3),in)     - updated coordinates
c     ks       (dp,sc,in)        - first face coordinate of projection
c     et       (dp,sc,in)        - second face coordinate of projection
c
c  output arguments:
c     tr       (dp,ar(3,3),out)  - transformation matrix for projection point
c
#include "impcom.inc"
c
      integer   nfacn,k
      double precision
     & x(4),y(4),z(4),ks,et,tr(3,3),ds(3),dt(3),vn(3),dfi(2,8),
     & ZERO,ONE,QUART,dxdksi,dxdeta,dydksi,dydeta,dzdksi,dzdeta,drdksi

      parameter (QUART=0.25d0, ZERO=0.0d0, ONE=1.0d0)

      dfi(1,1)= -QUART*(ONE-et)
      dfi(1,2)=  QUART*(ONE-et)
      dfi(1,3)=  QUART*(ONE+et)
      dfi(1,4)= -QUART*(ONE+et)
      dfi(2,1)= -QUART*(ONE-ks)
      dfi(2,2)= -QUART*(ONE+ks)
      dfi(2,3)=  QUART*(ONE+ks)
      dfi(2,4)=  QUART*(ONE-ks)

      dxdksi = ZERO
      dxdeta = ZERO
      dydksi = ZERO
      dydeta = ZERO
      dzdksi = ZERO
      dzdeta = ZERO
      do k = 1,nfacn
         dxdksi = dxdksi+dfi(1,k)*x(k)
         dxdeta = dxdeta+dfi(2,k)*x(k)
         dydksi = dydksi+dfi(1,k)*y(k)
         dydeta = dydeta+dfi(2,k)*y(k)
         dzdksi = dzdksi+dfi(1,k)*z(k)
         dzdeta = dzdeta+dfi(2,k)*z(k)
      enddo
      
      drdksi = dxdksi*dxdksi+dydksi*dydksi+dzdksi*dzdksi
      drdksi = sqrt(drdksi)
      ds(1) = dxdksi/drdksi
      ds(2) = dydksi/drdksi
      ds(3) = dzdksi/drdksi

c --- compute normal vector
      call v_norm (vn(1),3)
      call v_norm (ds(1),3)
      call v_cross (vn(1),ds(1),dt(1))
      call v_norm (dt(1),3)
c --- compute average tangent vector in general s direction (as t x n)
      call v_cross (dt(1),vn(1),ds(1))

c --- normalize vectors for transformation matrix
      call v_norm (ds(1),3)

c --- row-by-row fill of tr matrix
      tr(1,1) = ds(1)
      tr(1,2) = ds(2)
      tr(1,3) = ds(3)
      tr(2,1) = -dt(1)
      tr(2,2) = -dt(2)
      tr(2,3) = -dt(3)
      tr(3,1) = -vn(1)
      tr(3,2) = -vn(2)
      tr(3,3) = -vn(3)
c
 999  continue
      return
      end