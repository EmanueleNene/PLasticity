      subroutine InterpHygStrn(x,y,z,NdEps,Eps)
c --- primary function: Uses the nodal(corner) strain values to calculate strain at X (assuming linear shape functions)
c
c Input arguments:
c  x,y,z   (dp,sca,in)                - Point position in natural coordinates
c  NdEps   (dp,arr(8),in)             - Strains at nodes (corner points)

c
c output arguments:
c  Eps     (dp,sca,out)               - Strain at point (x,y,z)
c
#include "impcom.inc"
c


      integer             H_nElNode
      double precision    x,y,z, NdEps(8),Eps

      double precision    N(8)
      double precision    xi(3,8)
      
      double precision    v_dot

      integer             i

      double precision    ONE
      parameter           (ONE = 1.0d0)

      call NdOrder(xi)

      do i=1,8
        N(i)= (ONE+xi(1,i)*x)*(ONE+xi(2,i)*y)*(ONE+xi(3,i)*z)*0.125d0
      enddo

      Eps = v_dot(N(1),NdEps(1),8)


999   return
      end
