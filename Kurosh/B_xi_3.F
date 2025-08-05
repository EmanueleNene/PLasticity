      subroutine B_xi_3(x,y,z,H_nElNode,Bxi)
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
      double precision    x,y,z
      double precision    Bxi(3,H_nElNode)
      
      double precision    One, Two, Hlf, Qrt, Egth
      parameter          (One  = 1.000d0)
      parameter          (Two  = 2.000d0)
      parameter          (Qrt  = 0.250d0)
      parameter          (Hlf  = 0.500d0)
      parameter          (Egth = 0.125d0)

      if (H_nElNode .eq. 8) then
     
          Bxi(1,1) = -Egth * (One-y) * (One-z)
          Bxi(1,2) =  Egth * (One-y) * (One-z)
          Bxi(1,3) =  Egth * (One+y) * (One-z)
          Bxi(1,4) = -Egth * (One+y) * (One-z)
          Bxi(1,5) = -Egth * (One-y) * (One+z)
          Bxi(1,6) =  Egth * (One-y) * (One+z)
          Bxi(1,7) =  Egth * (One+y) * (One+z)
          Bxi(1,8) = -Egth * (One+y) * (One+z)
          
          Bxi(2,1) = -Egth * (One-x) * (One-z)
          Bxi(2,2) = -Egth * (One+x) * (One-z)
          Bxi(2,3) =  Egth * (One+x) * (One-z)
          Bxi(2,4) =  Egth * (One-x) * (One-z)
          Bxi(2,5) = -Egth * (One-x) * (One+z)
          Bxi(2,6) = -Egth * (One+x) * (One+z)
          Bxi(2,7) =  Egth * (One+x) * (One+z)
          Bxi(2,8) =  Egth * (One-x) * (One+z)
          
          Bxi(3,1) = -Egth * (One-x) * (One-y)
          Bxi(3,2) = -Egth * (One+x) * (One-y)
          Bxi(3,3) = -Egth * (One+x) * (One+y)
          Bxi(3,4) = -Egth * (One-x) * (One+y)
          Bxi(3,5) =  Egth * (One-x) * (One-y)
          Bxi(3,6) =  Egth * (One+x) * (One-y)
          Bxi(3,7) =  Egth * (One+x) * (One+y)
          Bxi(3,8) =  Egth * (One-x) * (One+y)

      elseif (H_nElNode .eq. 20) then
          Bxi(1,1 ) =  Egth* (y-1) * (z-1) * ( Two*x+y+z+1)
          Bxi(1,2 ) = -Egth* (y-1) * (z-1) * (-Two*x+y+z+1)
          Bxi(1,3 ) = -Egth* (y+1) * (z-1) * ( Two*x+y-z-1)
          Bxi(1,4 ) = -Egth* (y+1) * (z-1) * ( Two*x-y+z+1)
          Bxi(1,5 ) = -Egth* (y-1) * (z+1) * ( Two*x+y-z+1)
          Bxi(1,6 ) = -Egth* (y-1) * (z+1) * ( Two*x-y+z-1)
          Bxi(1,7 ) =  Egth* (y+1) * (z+1) * ( Two*x+y+z-1)
          Bxi(1,8 ) =  Egth* (y+1) * (z+1) * ( Two*x-y-z+1)
          Bxi(1,9 ) =  Qrt* (y*y-1) * (z-1) 
          Bxi(1,10) =  Hlf*x* (y+1) * (z-1) 
          Bxi(1,11) = -Qrt* (y*y-1) * (z-1) 
          Bxi(1,12) = -Hlf*x* (y-1) * (z-1) 
          Bxi(1,13) = -Qrt* (y*y-1) * (z+1) 
          Bxi(1,14) = -Hlf*x* (y+1) * (z+1) 
          Bxi(1,15) =  Qrt* (y*y-1) * (z+1) 
          Bxi(1,16) =  Hlf*x* (y-1) * (z+1) 
          Bxi(1,17) =  Qrt* (z*z-1) * (y-1) 
          Bxi(1,18) = -Qrt* (z*z-1) * (y+1) 
          Bxi(1,19) =  Qrt* (z*z-1) * (y+1) 
          Bxi(1,20) = -Qrt* (z*z-1) * (y-1)
      
          Bxi(2,1 ) =  Egth* (x-1) * (z-1) * ( x+Two*y+z+1)
          Bxi(2,2 ) = -Egth* (x+1) * (z-1) * (-x+Two*y+z+1)
          Bxi(2,3 ) = -Egth* (x+1) * (z-1) * ( x+Two*y-z-1)
          Bxi(2,4 ) = -Egth* (x-1) * (z-1) * ( x-Two*y+z+1)
          Bxi(2,5 ) = -Egth* (x-1) * (z+1) * ( x+Two*y-z+1)
          Bxi(2,6 ) = -Egth* (x+1) * (z+1) * ( x-Two*y+z-1)
          Bxi(2,7 ) =  Egth* (x+1) * (z+1) * ( x+Two*y+z-1)
          Bxi(2,8 ) =  Egth* (x-1) * (z+1) * ( x-Two*y-z+1)
          Bxi(2,9 ) =  Hlf*y* (x+1) * (z-1)
          Bxi(2,10) =  Qrt* (x*x-1) * (z-1)
          Bxi(2,11) = -Hlf*y* (x-1) * (z-1)
          Bxi(2,12) = -Qrt* (x*x-1) * (z-1)
          Bxi(2,13) = -Hlf*y* (x+1) * (z+1)
          Bxi(2,14) = -Qrt* (x*x-1) * (z+1)
          Bxi(2,15) =  Hlf*y* (x-1) * (z+1)
          Bxi(2,16) =  Qrt* (x*x-1) * (z+1)
          Bxi(2,17) =  Qrt* (z*z-1) * (x+1)
          Bxi(2,18) = -Qrt* (z*z-1) * (x+1)
          Bxi(2,19) =  Qrt* (z*z-1) * (x-1)
          Bxi(2,20) = -Qrt* (z*z-1) * (x-1)
      
          Bxi(3,1 ) =   Egth* (x-1) * (y-1) *( x+y+Two*z+1)
          Bxi(3,2 ) =  -Egth* (x+1) * (y-1) *(-x+y+Two*z+1)
          Bxi(3,3 ) =  -Egth* (x+1) * (y+1) *( x+y-Two*z-1)
          Bxi(3,4 ) =  -Egth* (x-1) * (y+1) *( x-y+Two*z+1)
          Bxi(3,5 ) =  -Egth* (x-1) * (y-1) *( x+y-Two*z+1)
          Bxi(3,6 ) =  -Egth* (x+1) * (y-1) *( x-y+Two*z-1)
          Bxi(3,7 ) =   Egth* (x+1) * (y+1) *( x+y+Two*z-1)
          Bxi(3,8 ) =   Egth* (x-1) * (y+1) *( x-y-Two*z+1)
          Bxi(3,9 ) =   Qrt* (y*y-1) * (x+1)
          Bxi(3,10) =   Qrt* (x*x-1) * (y+1)
          Bxi(3,11) =  -Qrt* (y*y-1) * (x-1)
          Bxi(3,12) =  -Qrt* (x*x-1) * (y-1)
          Bxi(3,13) =  -Qrt* (y*y-1) * (x+1)
          Bxi(3,14) =  -Qrt* (x*x-1) * (y+1)
          Bxi(3,15) =   Qrt* (y*y-1) * (x-1)
          Bxi(3,16) =   Qrt* (x*x-1) * (y-1)
          Bxi(3,17) =   Hlf*z* (x+1) * (y-1)
          Bxi(3,18) =  -Hlf*z* (x+1) * (y+1)
          Bxi(3,19) =   Hlf*z* (x-1) * (y+1)
          Bxi(3,20) =  -Hlf*z* (x-1) * (y-1)
      else
          call erhandler ('BndHygStrn>B_xi',5000,2,
     &                'Wrong number of nodes for 2D element',0.0d0,' ')
      endif
      
999   return
      end      