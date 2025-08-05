      subroutine B_Fcn_3(Bxi,Jinv,H_nElNode,B)
c --- primary function: Find B matix in terms of local x and y
c
c Input arguments:
c  Bxi     (dp,arr(2,H_nElNode),out)  - B matrix in local coordinate system
c  Jinv    (dp,arr(Dim,Dim),out)      - inverse Jacobian matrix
c  Dim     (int,sc,in)                - Dimension of the element (2D or 3D)
c  H_nElNode (int,sc,in)              - Number of nodes per element

c
c output arguments:
c  B     (dp,arr(3,2*H_nElNode),out) - B matrix in global coordinate system
c
#include "impcom.inc"
c

      
      integer             H_nElNode
      double precision    Bxi(3,H_nElNode),Jinv(3,3)
      double precision    B(6,3*H_nElNode)!,B1(16)
      
      double precision    Bx(3,H_nElNode)
      integer             i,loc
      
      call AxB_C(Jinv,Bxi,Bx,3,3,3,3,H_nElNode,3)
      call v_zero (B(1,1),18*H_nElNode)
      
      
      loc = 1
      do i=1,H_nElNode
          B(1,loc)   = Bx(1,i)
          B(5,loc)   = Bx(3,i)
          B(6,loc)   = Bx(2,i)
          
          B(2,loc+1) = Bx(2,i)
          B(4,loc+1) = Bx(3,i)
          B(6,loc+1) = Bx(1,i)
          
          B(3,loc+2) = Bx(3,i)
          B(4,loc+2) = Bx(2,i)
          B(5,loc+2) = Bx(1,i)
          
          loc = loc + 3
      enddo
            
999   return
      end      