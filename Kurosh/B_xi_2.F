      subroutine B_xi_2(x,y,H_nElNode,Bxi)
c --- primary function: Find B matix in terms of xi, eta, etc.
c
c Input arguments:
c   x      (dp,sc,in)        - coordinate in local isoparametric system
c   y      (dp,sc,in)        - coordinate in local isoparametric system
c   H_nElNode (int,sc,in)       - Number of nodes per element
c
c output arguments:
c  Bxi    (dp,arr(2,H_nElNode),out) -B matrix
c
#include "impcom.inc"
c

      
      integer             H_nElNode
      double precision    x,y
      double precision    Bxi(2,8)
      
      double precision    xx,yy
      double precision    Qrt, Hlf
      parameter          (Qrt = 0.25d0)
      parameter          (Hlf = 0.5d0)

      if (H_nElNode .eq. 4) then
          xx = x * Qrt
          yy = y * Qrt
      
          Bxi(1,1) =  yy-Qrt
          Bxi(1,2) = -yy+Qrt
          Bxi(1,3) =  yy+Qrt
          Bxi(1,4) = -yy-Qrt
      
          Bxi(2,1) =  xx-Qrt
          Bxi(2,2) = -xx-Qrt
          Bxi(2,3) =  xx+Qrt
          Bxi(2,4) = -xx+Qrt
      elseif (H_nElNode .eq. 8) then

          Bxi(1,1) =  Qrt * (1-y) * (2*x + y)
          Bxi(1,2) =  Qrt * (1-y) * (2*x - y)
          Bxi(1,3) =  Qrt * (1+y) * (2*x + y)
          Bxi(1,4) =  Qrt * (1+y) * (2*x - y)
          Bxi(1,5) = -x  * (1-y)
          Bxi(1,6) =  Hlf * (1-y*y)
          Bxi(1,7) = -x  * (1+y)
          Bxi(1,8) = -Hlf * (1-y*y)
              
          Bxi(2,1) =  Qrt * (1-x) * (2*y + x)
          Bxi(2,2) =  Qrt * (1+x) * (2*y - x)
          Bxi(2,3) =  Qrt * (1+x) * (2*y + x)
          Bxi(2,4) =  Qrt * (1-x) * (2*y - x)
          Bxi(2,5) = -Hlf * (1-x*x)
          Bxi(2,6) = -y  * (1+x)
          Bxi(2,7) =  Hlf * (1-x*x)
          Bxi(2,8) = -y  * (1-x)
      else
          call erhandler ('BndHygStrn>B_xi',5000,2,
     &                'Wrong number of nodes for 2D element',0.0d0,' ')
      endif
      
999   return
      end      