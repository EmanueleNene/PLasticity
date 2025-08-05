      subroutine GetEqCoef (xa,ya,za,xb,yb,zb,alpha,locu,locv,
     &                      al,am,ar,bl,bm,br)
c --- primary function: solve system of quadratic equations
c                       coef(1)*x^2 + coef(2)*x + coef(3) = 0.
c
c     According to the magnitudes of the 3 coef and also the determinent
c     the number of roots are returned. if a coef is relatively small
c     comparing to the maximum coef, then we make the small coef 0.
c
c  input arguments:
c     (dp, coef(3),in) -  the coeficiects of (ax^2 + bx + c = 0).
c                         Here a = coef(1).
#include "impcom.inc"
c
      double precision xa(7),ya(7),za(7),xb(7),yb(7),zb(7),alpha,
     &                 La,Lb
      double precision al(4,3),am(4,3),ar(4,3),bl(4,3),bm(4,3),br(4,3)
      double precision a31,a32,a33,a34,a36,a21,a22,a23,a24,a26,a11,a13,
     &                 a14,a16,ZERO,HALF
      parameter (ZERO=0.d0,HALF=0.5d0)
      integer locu,locv
      data a34/-16.d0/,a36/42.d0/,a31/-22.d0/,a33/-18.d0/,a32/14.d0/
      data a24/-13.8d0/,a26/39.6d0/,a21/-44.d0/,a23/8.4d0/,a22/9.8d0/
      data a14/-3.6d0/,a16/10.8d0/,a11/-10.8d0/,a13/19.6d0/


      La = (xa(1)-xa(2))**2+(ya(1)-ya(2))**2+(za(1)-za(2))**2
      Lb = (xb(2)-xb(3))**2+(yb(2)-yb(3))**2+(zb(2)-zb(3))**2
      if (La.lt.Lb) La = Lb
      La = sqrt(La)
      
      am(1,1) = xa(3)/La
      am(2,1) = HALF*(xa(2)-xa(1))/La
      am(3,1) = HALF*(xa(1)+xa(2)-2.d0*xa(3))/La
      am(4,1) = ZERO
      am(1,2) = ya(3)/La
      am(2,2) = HALF*(ya(2)-ya(1))/La
      am(3,2) = HALF*(ya(1)+ya(2)-2.d0*ya(3))/La      
      am(4,2) = ZERO
      am(1,3) = za(3)/La
      am(2,3) = HALF*(za(2)-za(1))/La
      am(3,3) = HALF*(za(1)+za(2)-2.d0*za(3))/La
      am(4,3) = ZERO

      bm(1,1) = xb(3)/La
      bm(2,1) = HALF*(xb(2)-xb(1))/La
      bm(3,1) = HALF*(xb(1)+xb(2)-2.d0*xb(3))/La
      bm(4,1) = ZERO
      bm(1,2) = yb(3)/La
      bm(2,2) = HALF*(yb(2)-yb(1))/La
      bm(3,2) = HALF*(yb(1)+yb(2)-2.d0*yb(3))/La
      bm(4,2) = ZERO
      bm(1,3) = zb(3)/La
      bm(2,3) = HALF*(zb(2)-zb(1))/La
      bm(3,3) = HALF*(zb(1)+zb(2)-2.d0*zb(3))/La
      bm(4,3) = ZERO 
      
      La =16.d0*La

      if (locu.ne.-1.and.locu.ne.2) then
          al(4,1) = (-xa(4)+2.d0*(xa(6)-xa(3))+xa(2))/La/alpha
          al(3,1) = (a34*xa(4)+a36*xa(6)+a31*xa(1)+a33*xa(3)+a32*xa(2))/
     &               La

          al(2,1) = (a24*xa(4)+a26*xa(6)+a21*xa(1)+a23*xa(3)+a22*xa(2))/
     &               La

          al(1,1) = (a14*xa(4)+a16*xa(6)+a11*xa(1)+a13*xa(3))/La

          al(4,2) = (-ya(4)+2.d0*(ya(6)-ya(3))+ya(2))/La/alpha
          al(3,2) = (a34*ya(4)+a36*ya(6)+a31*ya(1)+a33*ya(3)+a32*ya(2))/
     &               La

          al(2,2) = (a24*ya(4)+a26*ya(6)+a21*ya(1)+a23*ya(3)+a22*ya(2))/
     &               La

          al(1,2) = (a14*ya(4)+a16*ya(6)+a11*ya(1)+a13*ya(3))/La

          al(4,3) = (-za(4)+2.d0*(za(6)-za(3))+za(2))/La/alpha
          al(3,3) = (a34*za(4)+a36*za(6)+a31*za(1)+a33*za(3)+a32*za(2))/
     &               La

          al(2,3) = (a24*za(4)+a26*za(6)+a21*za(1)+a23*za(3)+a22*za(2))/
     &               La

          al(1,3) = (a14*za(4)+a16*za(6)+a11*za(1)+a13*za(3))/La

      else
          call v_move (am(1,1),al(1,1),12)
      endif
      if (locu.ne.1.and.locu.ne.2) then
          ar(4,1) = (-xa(1)+2.d0*(xa(3)-xa(7))+xa(5))/La/alpha
          ar(3,1) = (a34*xa(5)+a36*xa(7)+a31*xa(2)+a33*xa(3)+a32*xa(1))/
     &               La

          ar(2,1) = (-a24*xa(5)-a26*xa(7)-a21*xa(2)-a23*xa(3)-a22*xa(1))
     &               /La

          ar(1,1) = (a14*xa(5)+a16*xa(7)+a11*xa(2)+a13*xa(3))/La

          ar(4,2) = (-ya(1)+2.d0*(ya(3)-ya(7))+ya(5))/La/alpha
          ar(3,2) = (a34*ya(5)+a36*ya(7)+a31*ya(2)+a33*ya(3)+a32*ya(1))/
     &              La

          ar(2,2) = (-a24*ya(5)-a26*ya(7)-a21*ya(2)-a23*ya(3)-a22*ya(1))
     &              /La

          ar(1,2) = (a14*ya(5)+a16*ya(7)+a11*ya(2)+a13*ya(3))/La

          ar(4,3) = (-za(1)+2.d0*(za(3)-za(7))+za(5))/La/alpha
          ar(3,3) = (a34*za(5)+a36*za(7)+a31*za(2)+a33*za(3)+a32*za(1))/
     &              La

          ar(2,3) = (-a24*za(5)-a26*za(7)-a21*za(2)-a23*za(3)-a22*za(1))
     &              /La

          ar(1,3) = (a14*za(5)+a16*za(7)+a11*za(2)+a13*za(3))/La
      else
          call v_move (am(1,1),ar(1,1),12)
      endif

      if (locv.ne.-1.and.locv.ne.2) then
          bl(4,1) = (-xb(4)+2.d0*(xb(6)-xb(3))+xb(2))/La/alpha
          bl(3,1) = (a34*xb(4)+a36*xb(6)+a31*xb(1)+a33*xb(3)+a32*xb(2))/
     &              La

          bl(2,1) = (a24*xb(4)+a26*xb(6)+a21*xb(1)+a23*xb(3)+a22*xb(2))/
     &               La

          bl(1,1) = (a14*xb(4)+a16*xb(6)+a11*xb(1)+a13*xb(3))/La

          bl(4,2) = (-yb(4)+2.d0*(yb(6)-yb(3))+yb(2))/La/alpha
          bl(3,2) = (a34*yb(4)+a36*yb(6)+a31*yb(1)+a33*yb(3)+a32*yb(2))/
     &               La

          bl(2,2) = (a24*yb(4)+a26*yb(6)+a21*yb(1)+a23*yb(3)+a22*yb(2))/
     &               La

          bl(1,2) = (a14*yb(4)+a16*yb(6)+a11*yb(1)+a13*yb(3))/La

          bl(4,3) = (-zb(4)+2.d0*(zb(6)-zb(3))+zb(2))/La/alpha
          bl(3,3) = (a34*zb(4)+a36*zb(6)+a31*zb(1)+a33*zb(3)+a32*zb(2))/
     &               La

          bl(2,3) = (a24*zb(4)+a26*zb(6)+a21*zb(1)+a23*zb(3)+a22*zb(2))/
     &               La

          bl(1,3) = (a14*zb(4)+a16*zb(6)+a11*zb(1)+a13*zb(3))/La

      else
          call v_move (bm(1,1),bl(1,1),12)
      endif
      if (locv.ne.1.and.locv.ne.2) then
          br(4,1) = (-xb(1)+2.d0*(xb(3)-xb(7))+xb(5))/La/alpha
          br(3,1) = (a34*xb(5)+a36*xb(7)+a31*xb(2)+a33*xb(3)+a32*xb(1))/
     &               La

          br(2,1) = (-a24*xb(5)-a26*xb(7)-a21*xb(2)-a23*xb(3)-a22*xb(1))
     &              /La

          br(1,1) = (a14*xb(5)+a16*xb(7)+a11*xb(2)+a13*xb(3))/La

          br(4,2) = (-yb(1)+2.d0*(yb(3)-yb(7))+yb(5))/La/alpha
          br(3,2) = (a34*yb(5)+a36*yb(7)+a31*yb(2)+a33*yb(3)+a32*yb(1))/
     &               La

          br(2,2) = (-a24*yb(5)-a26*yb(7)-a21*yb(2)-a23*yb(3)-a22*yb(1))
     &              /La

          br(1,2) = (a14*yb(5)+a16*yb(7)+a11*yb(2)+a13*yb(3))/La

          br(4,3) = (-zb(1)+2.d0*(zb(3)-zb(7))+zb(5))/La/alpha
          br(3,3) = (a34*zb(5)+a36*zb(7)+a31*zb(2)+a33*zb(3)+a32*zb(1))/
     &               La

          br(2,3) = (-a24*zb(5)-a26*zb(7)-a21*zb(2)-a23*zb(3)-a22*zb(1))
     &              /La

          br(1,3) = (a14*zb(5)+a16*zb(7)+a11*zb(2)+a13*zb(3))/La
      else
          call v_move (bm(1,1),br(1,1),12)
      endif
c
      return
      end
