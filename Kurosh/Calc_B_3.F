      subroutine Calc_B_3(XYc,H,B,W)
c     &                      H_nElNode)
c --- primary function: Calculate stiffness matrix and force vector
c
c Input arguments:
c   XYc    (dp,arr(4,2),in)  - Corner node corrdinates
c   EL     (dp,sc,in)  -Longitudianl elasticity modulus
c   ET     (dp,sc,in)  -Transverse elasticity modulus
c   G      (dp,sc,in)  -In-plane shear modulus (GLT)
c   NLT    (dp,sc,in)  -In-plane Poisson's ratio
c   NTL    (dp,sc,in)  -In-plane Poisson's ratio
c   H      (dp,sc,in)  -Height of cross-section 
c   BetaL  (dp,sc,in)  -Longitudianl hydroexpansion coefficient
c   BetaT  (dp,sc,in)  -Transverse hydroexpansion coefficient
c   Dim    (int,sc,in) -Problem dimension (2D/3D)
c   H_nElNode (int,sc,in) -Number of nodes per element
c   H_nGptx  (int,sc,in) -Number of Gauss points to be used for integration in x-direction
c   H_nGpty  (int,sc,in) -Number of Gauss points to be used for integration in y-direction
c
c output arguments:
c  K     (dp,arr(2*H_nElNode,2*H_nElNode),out) - stiffness matrix
c  F     (dp,arr(2*H_nElNode),out)       - force vector
c
      
#include "impcom.inc"
#include "beam_hyg.inc"


      double precision    XYc(4,2),H
      
      double precision    J(3,3),Jinv(3,3),detJ
      double precision    Bxi(3,H_nElNode)
      double precision    B(6,3*H_nElNode,H_nGpt**3)
      double precision    W(H_nGpt**3)
      
      double precision    XY(H_nElNode,3)
      
      integer             i,p,q, count

      double precision    Hlf
      parameter           (Hlf = 0.5d0)
      
      
      XY(1:4,1:2) = XYc
      XY(1:4,3)   = 0.0d0
      XY(5:8,1:2) = XYc
      XY(5:8,3)   =  H
      if (H_key(2) .eq. 2) then
          XY(9 ,1:2) = Hlf * (XYc(2,:)+XYc(3,:))
          XY(10,1:2) = Hlf * (XYc(3,:)+XYc(4,:))
          XY(11,1:2) = Hlf * (XYc(4,:)+XYc(1,:))
          XY(12,1:2) = Hlf * (XYc(1,:)+XYc(2,:))
          XY(9:12,3) = 0.0d0
          XY(13:16,1:2) = XY(9:12,1:2)
          XY(13:16,3)   =  H
          XY(17:19,1:2) = XYc(2:4,:)
          XY(20,1:2)    = XYc(1,:)
          XY(17:20,3)   = Hlf * H
      endif
            
      
      count = 0
      do i=1,H_nGpt
        do p=1,H_nGpt
          do q=1,H_nGpt
             
             count = count+1
             
             call B_xi_3(H_xG(i),H_xG(p),H_xG(q),H_nElNode,Bxi(1,1))
             call JacH(Bxi(1,1),XY(1,1),3,H_nElNode,
     &                J(1,1),Jinv(1,1),detJ)
             call B_Fcn_3(Bxi(1,1),Jinv(1,1),H_nElNode,B(1,1,count))
              
             W(count) = detJ*H_wG(i)*H_wG(p)*H_wG(q)
             
          enddo
        enddo  
      enddo

      return
      end