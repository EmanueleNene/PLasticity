      function ContCheck4 (nfacn,x0,y0,nx,ny,nz,x,y,z,
     &                     dmax,xs,ys,zs,dd,ipass)
c --- primary function: check for false contact
c
#include "impcom.inc"
c
      double precision x0,y0,nx(8),ny(8),nz(8),x(8),y(8),z(8),
     &                 e01,e02,e03,ContCheck4,dd,xs,ys,zs,dmax,norm,
     &                 ZERO,ONE,QUART,fi(8),xx,yy,zz,lb
      integer nfacn,k,ipass
c
      parameter (ZERO=0.0d0, ONE=1.0d0, QUART=0.25d0)
c
      dd = 1.01d0*dmax
      lb = 10.0d0
      if (ipass.le.-1) lb = 900
      if (abs(x0).gt.lb.or.abs(y0).gt.lb) goto 999

      fi(1)= QUART*(ONE-x0)*(ONE-y0)
      fi(2)= QUART*(ONE+x0)*(ONE-y0)
      fi(3)= QUART*(ONE+x0)*(ONE+y0)
      fi(4)= QUART*(ONE-x0)*(ONE+y0)

      e01 = ZERO
      e02 = ZERO
      e03 = ZERO
      xx = -xs
      yy = -ys
      zz = -zs
      do k = 1,nfacn
         xx = xx + fi(k)*x(k)
         yy = yy + fi(k)*y(k)
         zz = zz + fi(k)*z(k)
         e01 = e01 + fi(k)*nx(k)
         e02 = e02 + fi(k)*ny(k)
         e03 = e03 + fi(k)*nz(k)
      enddo
      
      norm = sqrt(e01*e01+e02*e02+e03*e03)
      if (norm.lt.TINY) then
          dd = dmax
      else
          dd = -(e01*xx + e02*yy + e03*zz)/norm
      endif
        
 999  continue
      ContCheck4 = -ONE
      if (abs(dd).ge.dmax) ContCheck4 = ONE

      return
      end

