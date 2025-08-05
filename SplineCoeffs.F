      function SplineCoeffs(XYZ,n,C)
c --- primary function: calculate the Coefficients of splines passing through n points (X,Y,Z)
c     (the method is in fact using tailor's approximation and gauss elemination)
c     x= C(1) + C(2)* t + C(3)* t^2 + C(4)* t^3
c     y= C(5) + C(6)* t + C(7)* t^2 + C(8)* t^3
c     z= C(9) + C(10)*t + C(11)*t^2 + C(12)*t^3
c     where:  "0 =< t =< 1"   for each segment
c     and x,y,z define the shape of spline

c     X,Y,Z :   Point Coordinates
c     n     :   Number of Points
c     C     :   Output Coefficients

#include "impcom.inc"
#include "constants.inc"
c
      PTRFTN  MemAlloc
      
      integer   n, SplineCoeffs
      double precision  XYZ(3,n), C(12,n-1)
      
      double precision X(n), Y(n), Z(n)
      double precision a(n), bx(n),   by(n),  bz(n),
     &                       dx(n),   dy(n),  dz(n),
     &                       ex(n),   ey(n),  ez(n),
     &                       sx(n),   sy(n),  sz(n),
     &                      s2x(n),  s2y(n), s2z(n)
      
      integer           i,j
      double precision  h

c****************************************************
      SplineCoeffs = -1

      X=XYZ(1,:)
      Y=XYZ(2,:)
      Z=XYZ(3,:)
      
      call v_zero(C(1,1), (n-1)*12)
      if (n .le. 1) then
        goto 999
      elseif (n .eq. 2) then
        C(1 ,1) = X(1)
        C(5 ,1) = Y(1)
        C(9 ,1) = Z(1)
        
        C(2 ,1) = X(2) - X(1)
        C(6 ,1) = Y(2) - Y(1)
        C(10,1) = Z(2) - Z(1)
        
        SplineCoeffs = 0
        goto 999
      endif
      
      call v_zero(  a(1), n)
      call v_zero( bx(1), n)
      call v_zero( by(1), n)
      call v_zero( bz(1), n)
      call v_zero( dx(1), n)
      call v_zero( dy(1), n)
      call v_zero( dz(1), n)
      call v_zero( ex(1), n)
      call v_zero( ey(1), n)
      call v_zero( ez(1), n)
      call v_zero( sx(1), n)
      call v_zero( sy(1), n)
      call v_zero( sz(1), n)
      call v_zero(s2x(1), n)
      call v_zero(s2y(1), n)
      call v_zero(s2z(1), n)

c first calculating slopes at nodal points
      a(2) = 4.0d0
      DO j = 3,n
        a(j) = 4.0d0 - 1.0d0/a(j-1);
      ENDDO

      h = 1.0d0/(n-1)
      DO j = 2,n
        ex(j) = X(j)-X(j-1)
        ey(j) = Y(j)-Y(j-1)
        ez(j) = Z(j)-Z(j-1)
      ENDDO

      dx(2) = 6.0d0*(ex(3)-ex(2))/h**2-s2x(1);
      dy(2) = 6.0d0*(ey(3)-ey(2))/h**2-s2y(1);
      dz(2) = 6.0d0*(ez(3)-ez(2))/h**2-s2z(1);
      DO j=3,n-2
        dx(j) = 6.0d0*(ex(j+1)-ex(j))/h**2.0d0;
        dy(j) = 6.0d0*(ey(j+1)-ey(j))/h**2.0d0;
        dz(j) = 6.0d0*(ez(j+1)-ez(j))/h**2.0d0;
      ENDDO
      IF (n .GT. 2) THEN
       dx(n-1) = 6.0d0*(ex(n)-ex(n-1))/h**2.0d0-s2x(n);
       dy(n-1) = 6.0d0*(ey(n)-ey(n-1))/h**2.0d0-s2y(n);
       dz(n-1) = 6.0d0*(ez(n)-ez(n-1))/h**2.0d0-s2z(n);

       bx(2) = dx(2);
       by(2) = dy(2);
       bz(2) = dz(2);
       DO j = 3,n-1
         bx(j) = dx(j)-bx(j-1)/a(j-1);
         by(j) = dy(j)-by(j-1)/a(j-1);
         bz(j) = dz(j)-bz(j-1)/a(j-1);
       ENDDO

       s2x(n-1) = bx(n-1)/a(n-1);
       s2y(n-1) = by(n-1)/a(n-1);
       s2z(n-1) = bz(n-1)/a(n-1);
      ENDIF
      DO j = n-2,2,-1
        s2x(j) = (bx(j)-s2x(j+1))/a(j);
        s2y(j) = (by(j)-s2y(j+1))/a(j);
        s2z(j) = (bz(j)-s2z(j+1))/a(j);
      ENDDO

      sx(1) = ex(2)/h-s2x(1)*h/3.0d0-s2x(2)*h/6.0d0;
      sy(1) = ey(2)/h-s2y(1)*h/3.0d0-s2y(2)*h/6.0d0;
      sz(1) = ez(2)/h-s2z(1)*h/3.0d0-s2z(2)*h/6.0d0;
      DO j = 2,n
        sx(j) = sx(j-1) + (s2x(j-1)+s2x(j))*h/2.0d0;
        sy(j) = sy(j-1) + (s2y(j-1)+s2y(j))*h/2.0d0;
        sz(j) = sz(j-1) + (s2z(j-1)+s2z(j))*h/2.0d0;
      ENDDO
c **************************************************************
c  finding coefficients using calculated slopes!
      DO j = 1,n-1
      
        C(1 ,j) = X(j)
        C(5 ,j) = Y(j)
        C(9 ,j) = Z(j)
        
        C(2 ,j) = sx(j) * h
        C(6 ,j) = sy(j) * h
        C(10,j) = sz(j) * h
        
        C(3 ,j) = 3.0d0*( X(j+1)-X(j) ) - 2.0d0 *sx(j)*h - sx(j+1)*h
        C(7 ,j) = 3.0d0*( Y(j+1)-Y(j) ) - 2.0d0 *sy(j)*h - sy(j+1)*h
        C(11,j) = 3.0d0*( Z(j+1)-Z(j) ) - 2.0d0 *sz(j)*h - sz(j+1)*h
        
        C(4 ,j) = 2.0d0*( X(j)-X(j+1) ) + sx(j)*h + sx(j+1)*h
        C(8 ,j) = 2.0d0*( Y(j)-Y(j+1) ) + sy(j)*h + sy(j+1)*h
        C(12,j) = 2.0d0*( Z(j)-Z(j+1) ) + sz(j)*h + sz(j+1)*h

      ENDDO
      
      SplineCoeffs = 0
999   return        
      end