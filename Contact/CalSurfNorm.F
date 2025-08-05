      subroutine CalSurfNorm (ks,et,e1,e2,e3,x,y,z,ipass)
c --- primary function: calculate the normal at the contact point
c
#include "impcom.inc"
c
      double precision ks,et,e1,e2,e3,x(4),y(4),z(4),ZERO,ONE,QUART,
     &  dfi(2,4),dxdksi,dxdeta,dydksi,dydeta,dzdksi,dzdeta,dtdx,dtdy,
     &  dtdz,dajc,lb
      integer k,ipass
c
      parameter (QUART=0.25d0, ZERO=0.d0, ONE=1.d0)
c
      lb = 10.d0
      if (ipass.le.-1) lb = 900.d0
      if (abs(ks).gt.lb.or.abs(et).gt.lb) then
          e1 = ZERO
          e2 = ZERO
          e3 = ZERO
          goto 999
      endif
      dfi(1,1) = -QUART*(ONE-et)
      dfi(1,2) =  QUART*(ONE-et)
      dfi(1,3) =  QUART*(ONE+et)
      dfi(1,4) = -QUART*(ONE+et)
      dfi(2,1) = -QUART*(ONE-ks)
      dfi(2,2) = -QUART*(ONE+ks)
      dfi(2,3) =  QUART*(ONE+ks)
      dfi(2,4) =  QUART*(ONE-ks)
      dxdksi = ZERO
      dxdeta = ZERO
      dydksi = ZERO
      dydeta = ZERO
      dzdksi = ZERO
      dzdeta = ZERO
      do k = 1,4
         dxdksi = dxdksi+dfi(1,k)*x(k)
         dxdeta = dxdeta+dfi(2,k)*x(k)
         dydksi = dydksi+dfi(1,k)*y(k)
         dydeta = dydeta+dfi(2,k)*y(k)
         dzdksi = dzdksi+dfi(1,k)*z(k)
         dzdeta = dzdeta+dfi(2,k)*z(k)
      enddo
c
      dtdx = (dydksi*dzdeta -dydeta*dzdksi)
      dtdy = (dzdksi*dxdeta -dzdeta*dxdksi)
      dtdz = (dxdksi*dydeta -dxdeta*dydksi)
      dajc = dtdx*dtdx+dtdy*dtdy+dtdz*dtdz
      if (dajc.le.ZERO) then
          e1 = ZERO
          e2 = ZERO
          e3 = ZERO
      else
          dajc = sqrt(dajc)
          e1 = dtdx/dajc
          e2 = dtdy/dajc
          e3 = dtdz/dajc
      endif
 
 999  continue
 
      return
      end
