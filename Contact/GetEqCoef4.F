      subroutine GetEqCoef4(x,y,z,nx,ny,nz,off,ndir,itry,coef)
c --- primary function: solve coef(1)*x^2 + coef(2)*x + coef(3) = 0.
c
#include "impcom.inc"
c
      double precision x(4),y(4),z(4),nx(4),ny(4),nz(4),off(3),
     &  xsum(4),ysum(4),zsum(4),coef(16,2),
     &  snx(4),sny(4),snz(4),snxx(9),snyy(9),snzz(9),snxy(9),snxz(9),
     &  snyz(9),diameter,diameter2,TWO
      parameter (TWO=2.0d0)
     
      integer ndir,itry

      diameter = (x(1)-x(2))**2+(y(1)-y(2))**2+(z(1)-z(2))**2
      diameter2 = (x(2)-x(3))**2+(y(2)-y(3))**2+(z(2)-z(3))**2
      diameter = max(diameter,diameter2)
      diameter2 = (x(3)-x(4))**2+(y(3)-y(4))**2+(z(3)-z(4))**2
      diameter = max(diameter,diameter2)
      diameter2 = (x(1)-x(4))**2+(y(1)-y(4))**2+(z(1)-z(4))**2
      diameter = max(diameter,diameter2)
      
      diameter = sqrt(diameter)

      xsum(1) = (x(1)-x(2)+x(3)-x(4))/diameter 
      xsum(2) = (-x(1)+x(2)+x(3)-x(4))/diameter 
      xsum(3) = (-x(1)-x(2)+x(3)+x(4))/diameter 
      xsum(4) = (x(1)+x(2)+x(3)+x(4)-4.d0*off(1))/diameter

      ysum(1) = (y(1)-y(2)+y(3)-y(4))/diameter 
      ysum(2) = (-y(1)+y(2)+y(3)-y(4))/diameter 
      ysum(3) = (-y(1)-y(2)+y(3)+y(4))/diameter 
      ysum(4) = (y(1)+y(2)+y(3)+y(4)-4.d0*off(2))/diameter       

      zsum(1) = (z(1)-z(2)+z(3)-z(4))/diameter 
      zsum(2) = (-z(1)+z(2)+z(3)-z(4))/diameter 
      zsum(3) = (-z(1)-z(2)+z(3)+z(4))/diameter 
      zsum(4) = (z(1)+z(2)+z(3)+z(4)-4.d0*off(3))/diameter       

      snx(1) = nx(1) - nx(2) + nx(3) - nx(4)
      snx(2) = -nx(1) + nx(2) + nx(3) - nx(4)
      snx(3) = -nx(1) - nx(2) + nx(3) + nx(4)
      snx(4) = nx(1) + nx(2) + nx(3) + nx(4)

      sny(1) = ny(1) - ny(2) + ny(3) - ny(4)
      sny(2) = -ny(1) + ny(2) + ny(3) - ny(4)
      sny(3) = -ny(1) - ny(2) + ny(3) + ny(4)
      sny(4) = ny(1) + ny(2) + ny(3) + ny(4)

      snz(1) = nz(1) - nz(2) + nz(3) - nz(4)
      snz(2) = -nz(1) + nz(2) + nz(3) - nz(4)
      snz(3) = -nz(1) - nz(2) + nz(3) + nz(4)
      snz(4) = nz(1) + nz(2) + nz(3) + nz(4)

      if (ndir.eq.1) then
          coef(1,1) = ysum(1)*snz(1)-zsum(1)*sny(1)
          coef(2,1) = ysum(1)*snz(2)+ysum(2)*snz(1) -
     &                (zsum(1)*sny(2)+zsum(2)*sny(1))
          coef(3,1) = ysum(1)*snz(3)+ysum(3)*snz(1) -
     &                (zsum(1)*sny(3)+zsum(3)*sny(1))
          coef(4,1) = ysum(1)*snz(4)+ysum(2)*snz(3)+
     &                ysum(3)*snz(2)+ysum(4)*snz(1) -
     &                (zsum(1)*sny(4)+zsum(2)*sny(3)) -
     &                (zsum(3)*sny(2)+zsum(4)*sny(1))
          coef(5,1) = ysum(2)*snz(2)-zsum(2)*sny(2)
          coef(6,1) = ysum(3)*snz(3)-zsum(3)*sny(3)  
          coef(7,1) = ysum(2)*snz(4)+ysum(4)*snz(2)-
     &                (zsum(2)*sny(4)+zsum(4)*sny(2))
          coef(8,1) = ysum(3)*snz(4)+ysum(4)*snz(3)-
     &                (zsum(3)*sny(4)+zsum(4)*sny(3))
          coef(9,1) = ysum(4)*snz(4)-zsum(4)*sny(4)

      elseif (ndir.eq.2) then

          coef(1,1) = xsum(1)*snz(1)-zsum(1)*snx(1)
          coef(2,1) = xsum(1)*snz(2)+xsum(2)*snz(1) -
     &                (zsum(1)*snx(2)+zsum(2)*snx(1))
          coef(3,1) = xsum(1)*snz(3)+xsum(3)*snz(1) -
     &                (zsum(1)*snx(3)+zsum(3)*snx(1))
          coef(4,1) = xsum(1)*snz(4)+xsum(2)*snz(3)+
     &                xsum(3)*snz(2)+xsum(4)*snz(1) -
     &                (zsum(1)*snx(4)+zsum(2)*snx(3)) -
     &                (zsum(3)*snx(2)+zsum(4)*snx(1))
          coef(5,1) = xsum(2)*snz(2)-zsum(2)*snx(2)
          coef(6,1) = xsum(3)*snz(3)-zsum(3)*snx(3)  
          coef(7,1) = xsum(2)*snz(4)+xsum(4)*snz(2)-
     &               (zsum(2)*snx(4)+zsum(4)*snx(2))
          coef(8,1) = xsum(3)*snz(4)+xsum(4)*snz(3)-
     &               (zsum(3)*snx(4)+zsum(4)*snx(3))
          coef(9,1) = xsum(4)*snz(4)-zsum(4)*snx(4)

      elseif (ndir.eq.3) then

          coef(1,1) = xsum(1)*sny(1)-ysum(1)*snx(1)
          coef(2,1) = xsum(1)*sny(2)+xsum(2)*sny(1) -
     &               (ysum(1)*snx(2)+ysum(2)*snx(1))
          coef(3,1) = xsum(1)*sny(3)+xsum(3)*sny(1) -
     &               (ysum(1)*snx(3)+ysum(3)*snx(1))
          coef(4,1) = xsum(1)*sny(4)+xsum(2)*sny(3)+
     &                xsum(3)*sny(2)+xsum(4)*sny(1) -
     &               (ysum(1)*snx(4)+ysum(2)*snx(3)) -
     &               (ysum(3)*snx(2)+ysum(4)*snx(1))
          coef(5,1) = xsum(2)*sny(2)-ysum(2)*snx(2)
          coef(6,1) = xsum(3)*sny(3)-ysum(3)*snx(3)  
          coef(7,1) = xsum(2)*sny(4)+xsum(4)*sny(2)-
     &               (ysum(2)*snx(4)+ysum(4)*snx(2))
          coef(8,1) = xsum(3)*sny(4)+xsum(4)*sny(3)-
     &               (ysum(3)*snx(4)+ysum(4)*snx(3))
          coef(9,1) = xsum(4)*sny(4)-ysum(4)*snx(4)

      endif

      if (itry.eq.1) goto 100
      if (ndir.eq.1) then

          sny(1) = sny(1)+snz(1)
          sny(2) = sny(2)+snz(2)
          sny(3) = sny(3)+snz(3)
          sny(4) = sny(4)+snz(4)
          ysum(1) = ysum(1)+zsum(1)
          ysum(2) = ysum(2)+zsum(2)
          ysum(3) = ysum(3)+zsum(3)
          ysum(4) = ysum(4)+zsum(4)

          coef(1,2) = xsum(1)*sny(1)-ysum(1)*snx(1)
          coef(2,2) = xsum(1)*sny(2)+xsum(2)*sny(1) -
     &               (ysum(1)*snx(2)+ysum(2)*snx(1))
          coef(3,2) = xsum(1)*sny(3)+xsum(3)*sny(1) -
     &               (ysum(1)*snx(3)+ysum(3)*snx(1))
          coef(4,2) = xsum(1)*sny(4)+xsum(2)*sny(3)+
     &                xsum(3)*sny(2)+xsum(4)*sny(1) -
     &               (ysum(1)*snx(4)+ysum(2)*snx(3)) -
     &               (ysum(3)*snx(2)+ysum(4)*snx(1))
          coef(5,2) = xsum(2)*sny(2)-ysum(2)*snx(2)
          coef(6,2) = xsum(3)*sny(3)-ysum(3)*snx(3)  
          coef(7,2) = xsum(2)*sny(4)+xsum(4)*sny(2)-
     &               (ysum(2)*snx(4)+ysum(4)*snx(2))
          coef(8,2) = xsum(3)*sny(4)+xsum(4)*sny(3)-
     &               (ysum(3)*snx(4)+ysum(4)*snx(3))
          coef(9,2) = xsum(4)*sny(4)-ysum(4)*snx(4)

      elseif (ndir.eq.2) then

          snz(1) = snx(1)+snz(1)
          snz(2) = snx(2)+snz(2)
          snz(3) = snx(3)+snz(3)
          snz(4) = snx(4)+snz(4)
          zsum(1) = xsum(1)+zsum(1)
          zsum(2) = xsum(2)+zsum(2)
          zsum(3) = xsum(3)+zsum(3)
          zsum(4) = xsum(4)+zsum(4)

          coef(1,2) = ysum(1)*snz(1)-zsum(1)*sny(1)
          coef(2,2) = ysum(1)*snz(2)+ysum(2)*snz(1) -
     &               (zsum(1)*sny(2)+zsum(2)*sny(1))
          coef(3,2) = ysum(1)*snz(3)+ysum(3)*snz(1) -
     &               (zsum(1)*sny(3)+zsum(3)*sny(1))
          coef(4,2) = ysum(1)*snz(4)+ysum(2)*snz(3)+
     &                ysum(3)*snz(2)+ysum(4)*snz(1) -
     &               (zsum(1)*sny(4)+zsum(2)*sny(3)) -
     &               (zsum(3)*sny(2)+zsum(4)*sny(1))
          coef(5,2) = ysum(2)*snz(2)-zsum(2)*sny(2)
          coef(6,2) = ysum(3)*snz(3)-zsum(3)*sny(3)  
          coef(7,2) = ysum(2)*snz(4)+ysum(4)*snz(2)-
     &               (zsum(2)*sny(4)+zsum(4)*sny(2))
          coef(8,2) = ysum(3)*snz(4)+ysum(4)*snz(3)-
     &               (zsum(3)*sny(4)+zsum(4)*sny(3))
          coef(9,2) = ysum(4)*snz(4)-zsum(4)*sny(4)

      elseif (ndir.eq.3) then

          snx(1) = sny(1)+snx(1)
          snx(2) = sny(2)+snx(2)
          snx(3) = sny(3)+snx(3)
          snx(4) = sny(4)+snx(4)
          xsum(1) = ysum(1)+xsum(1)
          xsum(2) = ysum(2)+xsum(2)
          xsum(3) = ysum(3)+xsum(3)
          xsum(4) = ysum(4)+xsum(4)

          coef(1,2) = xsum(1)*snz(1)-zsum(1)*snx(1)
          coef(2,2) = xsum(1)*snz(2)+xsum(2)*snz(1) -
     &               (zsum(1)*snx(2)+zsum(2)*snx(1))
          coef(3,2) = xsum(1)*snz(3)+xsum(3)*snz(1) -
     &               (zsum(1)*snx(3)+zsum(3)*snx(1))
          coef(4,2) = xsum(1)*snz(4)+xsum(2)*snz(3)+
     &                xsum(3)*snz(2)+xsum(4)*snz(1) -
     &               (zsum(1)*snx(4)+zsum(2)*snx(3)) -
     &               (zsum(3)*snx(2)+zsum(4)*snx(1))
          coef(5,2) = xsum(2)*snz(2)-zsum(2)*snx(2)
          coef(6,2) = xsum(3)*snz(3)-zsum(3)*snx(3)  
          coef(7,2) = xsum(2)*snz(4)+xsum(4)*snz(2)-
     &               (zsum(2)*snx(4)+zsum(4)*snx(2))
          coef(8,2) = xsum(3)*snz(4)+xsum(4)*snz(3)-
     &               (zsum(3)*snx(4)+zsum(4)*snx(3))
          coef(9,2) = xsum(4)*snz(4)-zsum(4)*snx(4)

      endif

      goto 900

 100  continue

      if (ndir.eq.1) then

          snxx(1) = -sny(1)*sny(1)-snz(1)*snz(1)
          snxx(2) = -TWO*(sny(1)*sny(2)+snz(1)*snz(2))
          snxx(3) = -TWO*(sny(1)*sny(3)+snz(1)*snz(3))
          snxx(4) = -sny(2)*sny(2)-snz(2)*snz(2)
          snxx(5) = -TWO*(sny(1)*sny(4)+sny(2)*sny(3)+
     &                    snz(1)*snz(4)+snz(2)*snz(3))
          snxx(6) = -sny(3)*sny(3)-snz(3)*snz(3)
          snxx(7) = -TWO*(sny(2)*sny(4)+snz(2)*snz(4))
          snxx(8) = -TWO*(sny(3)*sny(4)+snz(3)*snz(4))
          snxx(9) = -sny(4)*sny(4)-snz(4)*snz(4)

          snxy(1) = snx(1)*sny(1)
          snxy(2) = snx(1)*sny(2)+snx(2)*sny(1) 
          snxy(3) = snx(1)*sny(3)+snx(3)*sny(1)
          snxy(4) = snx(2)*sny(2)
          snxy(5) = snx(1)*sny(4)+snx(2)*sny(3)+
     &              snx(3)*sny(2)+snx(4)*sny(1)
          snxy(6) = snx(3)*sny(3)   
          snxy(7) = snx(2)*sny(4)+snx(4)*sny(2)
          snxy(8) = snx(3)*sny(4)+snx(4)*sny(3)
          snxy(9) = snx(4)*sny(4)

          snxz(1) = snx(1)*snz(1)
          snxz(2) = snx(1)*snz(2)+snx(2)*snz(1) 
          snxz(3) = snx(1)*snz(3)+snx(3)*snz(1)
          snxz(4) = snx(2)*snz(2)
          snxz(5) = snx(1)*snz(4)+snx(2)*snz(3)+
     &              snx(3)*snz(2)+snx(4)*snz(1)
          snxz(6) = snx(3)*snz(3)   
          snxz(7) = snx(2)*snz(4)+snx(4)*snz(2)
          snxz(8) = snx(3)*snz(4)+snx(4)*snz(3)
          snxz(9) = snx(4)*snz(4)

          coef(1,2) = snxx(1)*xsum(1)+snxy(1)*ysum(1)+snxz(1)*zsum(1)
          coef(2,2) = snxx(2)*xsum(1)+snxx(1)*xsum(2) + 
     &                snxy(2)*ysum(1)+snxy(1)*ysum(2) +
     &                snxz(2)*zsum(1)+snxz(1)*zsum(2)
          coef(3,2) = snxx(3)*xsum(1)+snxx(1)*xsum(3) + 
     &                snxy(3)*ysum(1)+snxy(1)*ysum(3) +
     &                snxz(3)*zsum(1)+snxz(1)*zsum(3)
          coef(4,2) = snxx(4)*xsum(1)+snxx(2)*xsum(2) + 
     &                snxy(4)*ysum(1)+snxy(2)*ysum(2) +
     &                snxz(4)*zsum(1)+snxz(2)*zsum(2)
          coef(5,2) = snxx(5)*xsum(1)+snxx(3)*xsum(2) + 
     &                snxy(5)*ysum(1)+snxy(3)*ysum(2) +
     &                snxz(5)*zsum(1)+snxz(3)*zsum(2) +
     &                snxx(2)*xsum(3)+snxx(1)*xsum(4) + 
     &                snxy(2)*ysum(3)+snxy(1)*ysum(4) + 
     &                snxz(2)*zsum(3)+snxz(1)*zsum(4)
          coef(6,2) = snxx(6)*xsum(1)+snxx(3)*xsum(3) + 
     &                snxy(6)*ysum(1)+snxy(3)*ysum(3) +
     &                snxz(6)*zsum(1)+snxz(3)*zsum(3)
          coef(7,2) = snxx(7)*xsum(1)+snxx(5)*xsum(2) + 
     &                snxy(7)*ysum(1)+snxy(5)*ysum(2) +
     &                snxz(7)*zsum(1)+snxz(5)*zsum(2) +
     &                snxx(4)*xsum(3)+snxx(2)*xsum(4) + 
     &                snxy(4)*ysum(3)+snxy(2)*ysum(4) + 
     &                snxz(4)*zsum(3)+snxz(2)*zsum(4)
          coef(8,2) = snxx(8)*xsum(1)+snxx(6)*xsum(2) + 
     &                snxy(8)*ysum(1)+snxy(6)*ysum(2) +
     &                snxz(8)*zsum(1)+snxz(6)*zsum(2) +
     &                snxx(5)*xsum(3)+snxx(3)*xsum(4) + 
     &                snxy(5)*ysum(3)+snxy(3)*ysum(4) + 
     &                snxz(5)*zsum(3)+snxz(3)*zsum(4)
          coef(9,2) = snxx(4)*xsum(2)+snxy(4)*ysum(2)+snxz(4)*zsum(2)
          coef(10,2) = snxx(6)*xsum(3)+snxy(6)*ysum(3)+snxz(6)*zsum(3)
          coef(11,2) = snxx(9)*xsum(1)+snxx(8)*xsum(2) + 
     &                 snxy(9)*ysum(1)+snxy(8)*ysum(2) +
     &                 snxz(9)*zsum(1)+snxz(8)*zsum(2) +
     &                 snxx(7)*xsum(3)+snxx(5)*xsum(4) + 
     &                 snxy(7)*ysum(3)+snxy(5)*ysum(4) + 
     &                 snxz(7)*zsum(3)+snxz(5)*zsum(4)
          coef(12,2) = snxx(7)*xsum(2)+snxx(4)*xsum(4) + 
     &                 snxy(7)*ysum(2)+snxy(4)*ysum(4) +
     &                 snxz(7)*zsum(2)+snxz(4)*zsum(4)
          coef(13,2) = snxx(8)*xsum(3)+snxx(6)*xsum(4) + 
     &                 snxy(8)*ysum(3)+snxy(6)*ysum(4) +
     &                 snxz(8)*zsum(3)+snxz(6)*zsum(4)
          coef(14,2) = snxx(9)*xsum(2)+snxx(7)*xsum(4) + 
     &                 snxy(9)*ysum(2)+snxy(7)*ysum(4) +
     &                 snxz(9)*zsum(2)+snxz(7)*zsum(4)
          coef(15,2) = snxx(9)*xsum(3)+snxx(8)*xsum(4) + 
     &                 snxy(9)*ysum(3)+snxy(8)*ysum(4) +
     &                 snxz(9)*zsum(3)+snxz(8)*zsum(4)
          coef(16,2) = snxx(9)*xsum(4)+snxy(9)*ysum(4)+snxz(9)*zsum(4)

      elseif (ndir.eq.2) then

          snxy(1) = snx(1)*sny(1)
          snxy(2) = snx(1)*sny(2)+snx(2)*sny(1) 
          snxy(3) = snx(1)*sny(3)+snx(3)*sny(1)
          snxy(4) = snx(2)*sny(2)
          snxy(5) = snx(1)*sny(4)+snx(2)*sny(3)+
     &              snx(3)*sny(2)+snx(4)*sny(1)
          snxy(6) = snx(3)*sny(3)   
          snxy(7) = snx(2)*sny(4)+snx(4)*sny(2)
          snxy(8) = snx(3)*sny(4)+snx(4)*sny(3)
          snxy(9) = snx(4)*sny(4)

          snyz(1) = snz(1)*sny(1)
          snyz(2) = snz(1)*sny(2)+snz(2)*sny(1) 
          snyz(3) = snz(1)*sny(3)+snz(3)*sny(1)
          snyz(4) = snz(2)*sny(2)
          snyz(5) = snz(1)*sny(4)+snz(2)*sny(3)+
     &              snz(3)*sny(2)+snz(4)*sny(1)
          snyz(6) = snz(3)*sny(3)   
          snyz(7) = snz(2)*sny(4)+snz(4)*sny(2)
          snyz(8) = snz(3)*sny(4)+snz(4)*sny(3)
          snyz(9) = snz(4)*sny(4)

          snyy(1) = -snx(1)*snx(1)-snz(1)*snz(1)
          snyy(2) = -TWO*(snx(1)*snx(2)+snz(1)*snz(2))
          snyy(3) = -TWO*(snx(1)*snx(3)+snz(1)*snz(3))
          snyy(4) = -snx(2)*snx(2)-snz(2)*snz(2)
          snyy(5) = -TWO*(snx(1)*snx(4)+snx(2)*snx(3)+
     &                    snz(1)*snz(4)+snz(2)*snz(3))
          snyy(6) = -snx(3)*snx(3)-snz(3)*snz(3)
          snyy(7) = -TWO*(snx(2)*snx(4)+snz(2)*snz(4))
          snyy(8) = -TWO*(snx(3)*snx(4)+snz(3)*snz(4))
          snyy(9) = -snx(4)*snx(4)-snz(4)*snz(4)

          coef(1,2) = snxy(1)*xsum(1)+snyy(1)*ysum(1)+snyz(1)*zsum(1)
          coef(2,2) = snxy(2)*xsum(1)+snxy(1)*xsum(2) + 
     &                snyy(2)*ysum(1)+snyy(1)*ysum(2) +
     &                snyz(2)*zsum(1)+snyz(1)*zsum(2)
          coef(3,2) = snxy(3)*xsum(1)+snxy(1)*xsum(3) + 
     &                snyy(3)*ysum(1)+snyy(1)*ysum(3) +
     &                snyz(3)*zsum(1)+snyz(1)*zsum(3)
          coef(4,2) = snxy(4)*xsum(1)+snxy(2)*xsum(2) + 
     &                snyy(4)*ysum(1)+snyy(2)*ysum(2) +
     &                snyz(4)*zsum(1)+snyz(2)*zsum(2)
          coef(5,2) = snxy(5)*xsum(1)+snxy(3)*xsum(2) + 
     &                snyy(5)*ysum(1)+snyy(3)*ysum(2) +
     &                snyz(5)*zsum(1)+snyz(3)*zsum(2) +
     &                snxy(2)*xsum(3)+snxy(1)*xsum(4) + 
     &                snyy(2)*ysum(3)+snyy(1)*ysum(4) + 
     &                snyz(2)*zsum(3)+snyz(1)*zsum(4)
          coef(6,2) = snxy(6)*xsum(1)+snxy(3)*xsum(3) + 
     &                snyy(6)*ysum(1)+snyy(3)*ysum(3) +
     &                snyz(6)*zsum(1)+snyz(3)*zsum(3)
          coef(7,2) = snxy(7)*xsum(1)+snxy(5)*xsum(2) + 
     &                snyy(7)*ysum(1)+snyy(5)*ysum(2) +
     &                snyz(7)*zsum(1)+snyz(5)*zsum(2) +
     &                snxy(4)*xsum(3)+snxy(2)*xsum(4) + 
     &                snyy(4)*ysum(3)+snyy(2)*ysum(4) + 
     &                snyz(4)*zsum(3)+snyz(2)*zsum(4)
          coef(8,2) = snxy(8)*xsum(1)+snxy(6)*xsum(2) + 
     &                snyy(8)*ysum(1)+snyy(6)*ysum(2) +
     &                snyz(8)*zsum(1)+snyz(6)*zsum(2) +
     &                snxy(5)*xsum(3)+snxy(3)*xsum(4) + 
     &                snyy(5)*ysum(3)+snyy(3)*ysum(4) + 
     &                snyz(5)*zsum(3)+snyz(3)*zsum(4)
          coef(9,2) = snxy(4)*xsum(2)+snyy(4)*ysum(2)+snyz(4)*zsum(2)
          coef(10,2) = snxy(6)*xsum(3)+snyy(6)*ysum(3)+snyz(6)*zsum(3)
          coef(11,2) = snxy(9)*xsum(1)+snxy(8)*xsum(2) + 
     &                 snyy(9)*ysum(1)+snyy(8)*ysum(2) +
     &                 snyz(9)*zsum(1)+snyz(8)*zsum(2) +
     &                 snxy(7)*xsum(3)+snxy(5)*xsum(4) + 
     &                 snyy(7)*ysum(3)+snyy(5)*ysum(4) + 
     &                 snyz(7)*zsum(3)+snyz(5)*zsum(4)
          coef(12,2) = snxy(7)*xsum(2)+snxy(4)*xsum(4) + 
     &                 snyy(7)*ysum(2)+snyy(4)*ysum(4) +
     &                 snyz(7)*zsum(2)+snyz(4)*zsum(4)
          coef(13,2) = snxy(8)*xsum(3)+snxy(6)*xsum(4) + 
     &                 snyy(8)*ysum(3)+snyy(6)*ysum(4) +
     &                 snyz(8)*zsum(3)+snyz(6)*zsum(4)
          coef(14,2) = snxy(9)*xsum(2)+snxy(7)*xsum(4) + 
     &                 snyy(9)*ysum(2)+snyy(7)*ysum(4) +
     &                 snyz(9)*zsum(2)+snyz(7)*zsum(4)
          coef(15,2) = snxy(9)*xsum(3)+snxy(8)*xsum(4) + 
     &                 snyy(9)*ysum(3)+snyy(8)*ysum(4) +
     &                 snyz(9)*zsum(3)+snyz(8)*zsum(4)
          coef(16,2) = snxy(9)*xsum(4)+snyy(9)*ysum(4)+snyz(9)*zsum(4)

      elseif (ndir.eq.3) then

          snxz(1) = snx(1)*snz(1)
          snxz(2) = snx(1)*snz(2)+snx(2)*snz(1) 
          snxz(3) = snx(1)*snz(3)+snx(3)*snz(1)
          snxz(4) = snx(2)*snz(2)
          snxz(5) = snx(1)*snz(4)+snx(2)*snz(3)+
     &              snx(3)*snz(2)+snx(4)*snz(1)
          snxz(6) = snx(3)*snz(3)   
          snxz(7) = snx(2)*snz(4)+snx(4)*snz(2)
          snxz(8) = snx(3)*snz(4)+snx(4)*snz(3)
          snxz(9) = snx(4)*snz(4)

          snyz(1) = snz(1)*sny(1)
          snyz(2) = snz(1)*sny(2)+snz(2)*sny(1) 
          snyz(3) = snz(1)*sny(3)+snz(3)*sny(1)
          snyz(4) = snz(2)*sny(2)
          snyz(5) = snz(1)*sny(4)+snz(2)*sny(3)+
     &              snz(3)*sny(2)+snz(4)*sny(1)
          snyz(6) = snz(3)*sny(3)   
          snyz(7) = snz(2)*sny(4)+snz(4)*sny(2)
          snyz(8) = snz(3)*sny(4)+snz(4)*sny(3)
          snyz(9) = snz(4)*sny(4)

          snzz(1) = -snx(1)*snx(1)-sny(1)*sny(1)
          snzz(2) = -TWO*(snx(1)*snx(2)+sny(1)*sny(2))
          snzz(3) = -TWO*(snx(1)*snx(3)+sny(1)*sny(3))
          snzz(4) = -snx(2)*snx(2)-sny(2)*sny(2)
          snzz(5) = -TWO*(snx(1)*snx(4)+snx(2)*snx(3)+
     &                     sny(1)*sny(4)+sny(2)*sny(3))
          snzz(6) = -snx(3)*snx(3)-sny(3)*sny(3)
          snzz(7) = -TWO*(snx(2)*snx(4)+sny(2)*sny(4))
          snzz(8) = -TWO*(snx(3)*snx(4)+sny(3)*sny(4))
          snzz(9) = -snx(4)*snx(4)-sny(4)*sny(4)

          coef(1,2) = snxz(1)*xsum(1)+snyz(1)*ysum(1)+snzz(1)*zsum(1)
          coef(2,2) = snxz(2)*xsum(1)+snxz(1)*xsum(2) + 
     &                snyz(2)*ysum(1)+snyz(1)*ysum(2) +
     &                snzz(2)*zsum(1)+snzz(1)*zsum(2)
          coef(3,2) = snxz(3)*xsum(1)+snxz(1)*xsum(3) + 
     &                snyz(3)*ysum(1)+snyz(1)*ysum(3) +
     &                snzz(3)*zsum(1)+snzz(1)*zsum(3)
          coef(4,2) = snxz(4)*xsum(1)+snxz(2)*xsum(2) + 
     &                snyz(4)*ysum(1)+snyz(2)*ysum(2) +
     &                snzz(4)*zsum(1)+snzz(2)*zsum(2)
          coef(5,2) = snxz(5)*xsum(1)+snxz(3)*xsum(2) + 
     &                snyz(5)*ysum(1)+snyz(3)*ysum(2) +
     &                snzz(5)*zsum(1)+snzz(3)*zsum(2) +
     &                snxz(2)*xsum(3)+snxz(1)*xsum(4) + 
     &                snyz(2)*ysum(3)+snyz(1)*ysum(4) + 
     &                snzz(2)*zsum(3)+snzz(1)*zsum(4)
          coef(6,2) = snxz(6)*xsum(1)+snxz(3)*xsum(3) + 
     &                snyz(6)*ysum(1)+snyz(3)*ysum(3) +
     &                snzz(6)*zsum(1)+snzz(3)*zsum(3)
          coef(7,2) = snxz(7)*xsum(1)+snxz(5)*xsum(2) + 
     &                snyz(7)*ysum(1)+snyz(5)*ysum(2) +
     &                snzz(7)*zsum(1)+snzz(5)*zsum(2) +
     &                snxz(4)*xsum(3)+snxz(2)*xsum(4) + 
     &                snyz(4)*ysum(3)+snyz(2)*ysum(4) + 
     &                snzz(4)*zsum(3)+snzz(2)*zsum(4)
          coef(8,2) = snxz(8)*xsum(1)+snxz(6)*xsum(2) + 
     &                snyz(8)*ysum(1)+snyz(6)*ysum(2) +
     &                snzz(8)*zsum(1)+snzz(6)*zsum(2) +
     &                snxz(5)*xsum(3)+snxz(3)*xsum(4) + 
     &                snyz(5)*ysum(3)+snyz(3)*ysum(4) + 
     &                snzz(5)*zsum(3)+snzz(3)*zsum(4)
          coef(9,2) = snxz(4)*xsum(2)+snyz(4)*ysum(2)+snzz(4)*zsum(2)
          coef(10,2) = snxz(6)*xsum(3)+snyz(6)*ysum(3)+snzz(6)*zsum(3)
          coef(11,2) = snxz(9)*xsum(1)+snxz(8)*xsum(2) + 
     &                 snyz(9)*ysum(1)+snyz(8)*ysum(2) +
     &                 snzz(9)*zsum(1)+snzz(8)*zsum(2) +
     &                 snxz(7)*xsum(3)+snxz(5)*xsum(4) + 
     &                 snyz(7)*ysum(3)+snyz(5)*ysum(4) + 
     &                 snzz(7)*zsum(3)+snzz(5)*zsum(4)
          coef(12,2) = snxz(7)*xsum(2)+snxz(4)*xsum(4) + 
     &                 snyz(7)*ysum(2)+snyz(4)*ysum(4) +
     &                 snzz(7)*zsum(2)+snzz(4)*zsum(4)
          coef(13,2) = snxz(8)*xsum(3)+snxz(6)*xsum(4) + 
     &                 snyz(8)*ysum(3)+snyz(6)*ysum(4) +
     &                 snzz(8)*zsum(3)+snzz(6)*zsum(4)
          coef(14,2) = snxz(9)*xsum(2)+snxz(7)*xsum(4) + 
     &                 snyz(9)*ysum(2)+snyz(7)*ysum(4) +
     &                 snzz(9)*zsum(2)+snzz(7)*zsum(4)
          coef(15,2) = snxz(9)*xsum(3)+snxz(8)*xsum(4) + 
     &                 snyz(9)*ysum(3)+snyz(8)*ysum(4) +
     &                 snzz(9)*zsum(3)+snzz(8)*zsum(4)
          coef(16,2) = snxz(9)*xsum(4)+snyz(9)*ysum(4)+snzz(9)*zsum(4)

      endif

 900  continue
      return
      end
