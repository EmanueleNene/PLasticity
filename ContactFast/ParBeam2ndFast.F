      subroutine ParBeam2ndFast (rad,xa,ya,za,xb,yb,zb,
     &                           dmax,epsil,nocon,tols)
c --- primary function: compute the relations between an integration point
c                       of the contact element and a 3d line segment
#include "impcom.inc"
#include "constants.inc"
c 
      integer          i,nocon
c
      double precision v_dot,rad(2),xa(3),ya(3),za(3),xb(3),yb(3),zb(3),
     &                 dmax,epsil,tols

      double precision ZERO,HALF,ONE,TWO,x0,y0,z0,
     &                 d,d0,d1,xmin,xmax,ymin,ymax,zmin,zmax,
     &                 a(3),ks,b(3),c,propen,e,cang,dang,angtol

c
      parameter (ZERO = 0.0d0, 
     &           HALF = 0.5d0,
     &           ONE = 1.0d0,
     &           TWO = 2.0d0,
     &           angtol = 0.d0)
c
      nocon = 2
      d = dmax
      
      a(1) = xa(3) - xa(1) 
      a(2) = ya(3) - ya(1) 
      a(3) = za(3) - za(1) 
      c = v_dot (a(1),a(1),3)
      
      if (angtol.gt.ZERO) then
          b(1) = xb(3) - xb(1) 
          b(2) = yb(3) - yb(1) 
          b(3) = zb(3) - zb(1) 
          e = v_dot (b(1),b(1),3)
          cang  = abs(v_dot(a(1),b(1),3)/sqrt(c)/sqrt(e))
          if (cang.lt.-ONE) cang = -ONE
          if (cang.gt.ONE) cang = ONE
          dang = acos(cang)*180.d0/W_Pi
          if (dang.gt.angtol) goto 998
      endif

c --- check out of pinball region
c 
      do i = 1,3,2
         x0 = xb(i)
         y0 = yb(i)
         z0 = zb(i)
      
         xmin = min(xa(1),xa(3)) - dmax
         if (x0.le.xmin) cycle
         xmax = max(xa(1),xa(3)) + dmax
         if (x0.ge.xmax) cycle
         ymin = min(ya(1),ya(3)) - dmax
         if (y0.le.ymin) cycle
         ymax = max(ya(1),ya(3)) + dmax
         if (y0.ge.ymax) cycle
         zmin = min(za(1),za(3)) - dmax
         if (z0.le.zmin) cycle
         zmax = max(za(1),za(3)) + dmax
         if (z0.ge.zmax) cycle

         b(1) = x0 - 0.5d0*(xa(3)+xa(1))
         b(2) = y0 - 0.5d0*(ya(3)+ya(1))
         b(3) = z0 - 0.5d0*(za(3)+za(1))
         d0 = TWO*v_dot (a(1),b(1),3)
         ks = d0/c
         propen = ONE + 0.01d0*tols
      
         if (abs(ks).gt.propen) cycle
      
         b(1) = (xa(1)-x0)
         b(2) = (ya(1)-y0)
         b(3) = (za(1)-z0)
         ks = HALF*(ks+ONE)
         d0 = (b(1)+a(1)*ks)**2+(b(2)+a(2)*ks)**2+(b(3)+a(3)*ks)**2
         d1 = abs(sqrt(d0) - abs(rad(1)) - abs(rad(2)))
         
         d = min(d,d1)
         
         nocon = 0
         
      enddo
 
998   continue  
      
      if (nocon.ne.0) d = dmax
      epsil = d
          
      return
      end

