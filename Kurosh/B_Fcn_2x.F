      subroutine B_Fcn_2x(Bxi,Jinv,H_nElNode,B1)
c --- primary function: Find B matix in x-direction in terms of local x and y
c
c Input arguments:
c  Bxi     (dp,arr(2,H_nElNode),out)  - B matrix in local coordinate system
c  Jinv    (dp,arr(Dim,Dim),out)      - inverse Jacobian matrix
c  Dim     (int,sc,in)                - Dimension of the element (2D or 3D)
c  H_nElNode (int,sc,in)              - Number of nodes per element

c
c output arguments:
c  B1     (dp,arr(3,2*H_nElNode),out) - B matrix in x-direction in global coordinate system
c
#include "impcom.inc"
c

      
      integer             H_nElNode
      double precision    Bxi(2,8),Jinv(2,2)
      double precision    B1(16)
      
      double precision    Bx(2,8)
      integer             i
      
      call AxB_C(Jinv,Bxi,Bx,2,2,2,2,H_nElNode,2)
      call v_zero (B1(1) ,2*H_nElNode)
      
      
      do i=1,H_nElNode
          B1(i*2-1)   = Bx(1,i)
      enddo
            
999   return
      end      