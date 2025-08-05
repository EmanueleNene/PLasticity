      subroutine JacH(Bxi,X,Dim,H_nElNode,J,Jinv,detJ)
c --- primary function: Find Jacobian Matrix
c
c Input arguments:
c  Bxi    (dp,arr(Dim,H_nElNode),out)  - B matrix in terms of isoparameters
c  X      (dp,arr(H_nElNode,Dim),in)   - Nodal corrdinates
c  Dim    (int,sc,in)               - Dimension of the element (2D or 3D)
c  H_nElNode (int,sc,in)               - Number of nodes per element
c
c output arguments:
c  J     (dp,arr(Dim,Dim),out) - Jacobian matrix
c  Jinv  (dp,arr(Dim,Dim),out) - inverse Jacobian matrix
c  detJ  (dp,sc,out)           - determinant of Jacobian matrix
c
#include "impcom.inc"
c

      
      integer             Dim, H_nElNode
      double precision    Bxi(Dim,H_nElNode),X(H_nElNode,Dim)
      double precision    J(Dim,Dim),Jinv(Dim,Dim),detJ
      
      double precision    dinv
      
      call AxB_C(Bxi,X,J,Dim,H_nElNode,Dim,Dim,Dim,H_nElNode)
      
      if (Dim .eq. 2) then
          detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1)
          dinv = 1.0d0/detJ
      
          Jinv(1,1) =  J(2,2)*dinv
          Jinv(2,2) =  J(1,1)*dinv
          Jinv(1,2) = -J(1,2)*dinv
          Jinv(2,1) = -J(2,1)*dinv
      elseif (Dim .eq. 3) then
          detJ = J(1,1) * ( J(2,2)*J(3,3) - J(2,3)*J(3,2) )
     &          -J(1,2) * ( J(2,1)*J(3,3) - J(2,3)*J(3,1) )
     &          +J(1,3) * ( J(2,1)*J(3,2) - J(2,2)*J(3,1) )
          
          dinv = 1.0d0/detJ
          
          Jinv(1,1) = dinv*( J(2,2)*J(3,3) - J(2,3)*J(3,2) )
          Jinv(1,2) = dinv*( J(1,3)*J(3,2) - J(1,2)*J(3,3) )
          Jinv(1,3) = dinv*( J(1,2)*J(2,3) - J(1,3)*J(2,2) )
          
          Jinv(2,1) = dinv*( J(2,3)*J(3,1) - J(2,1)*J(3,3) )
          Jinv(2,2) = dinv*( J(1,1)*J(3,3) - J(1,3)*J(3,1) )
          Jinv(2,3) = dinv*( J(1,3)*J(2,1) - J(1,1)*J(2,3) )
          
          Jinv(3,1) = dinv*( J(2,1)*J(3,2) - J(2,2)*J(3,1) )
          Jinv(3,2) = dinv*( J(1,2)*J(3,1) - J(1,1)*J(3,2) )
          Jinv(3,3) = dinv*( J(1,1)*J(2,2) - J(1,2)*J(2,1) )
      else
          call erhandler ('BndHygStrn>Jac',5000,2,
     &        'Wrong element dimension',0.0d0,' ')
      endif

999   return
      end      