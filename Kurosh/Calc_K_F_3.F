      subroutine Calc_K_F_3(XYc,EL,ET,G,NLT,NTL,H,BetaL,BetaT,DH,K,F)
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


      double precision    XYc(4,2),EL,ET,G,NLT,NTL,H,BetaL,BetaT,DH
      double precision    K(3*H_nElNode,3*H_nElNode),F(3*H_nElNode)

      double precision    J(3,3),Jinv(3,3),detJ
      double precision    Bxi(3,H_nElNode)
      double precision    B(6,3*H_nElNode)
      
      double precision    C(6,6),Beta(6)
      
      double precision    XY(H_nElNode,3)
      
      double precision    ff(3*H_nElNode,6),fff(3*H_nElNode,6)
      double precision    sss
      
      integer             i,p,q, count

      double precision    Hlf
      parameter           (Hlf = 0.5d0)
      
      Beta(1) = BetaL
      Beta(2) = BetaT
      Beta(3) = BetaT
      Beta(4) = 0.0d0
      Beta(5) = 0.0d0
      Beta(6) = 0.0d0
      
      call C_Mtx(EL,ET,G,G*0.1d0,NLT,NTL,0.0d0,6,C)
      
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
            
      call v_zero (K(1,1),9*H_nElNode*H_nElNode)
      call v_zero (F(1),3*H_nElNode)
      call v_zero (ff(1,1),18*H_nElNode)
      
      count = 0
      do i=1,H_nGpt
        do p=1,H_nGpt
          do q=1,H_nGpt
             call B_xi_3(H_xG(i),H_xG(p),H_xG(q),H_nElNode,Bxi(1,1))
             call JacH(Bxi(1,1),XY(1,1),3,H_nElNode,
     &                J(1,1),Jinv(1,1),detJ)
             call B_Fcn_3(Bxi(1,1),Jinv(1,1),H_nElNode,B(1,1))
              
             sss = detJ*H_wG(i)*H_wG(p)*H_wG(q)
             
             call cAtxBxA_C(B(1,1),C(1,1),K(1,1),6,6,
     &                 3*H_nElNode,6,3*H_nElNode,fff(1,1),sss)
             
             
             call v1pv2_v1(ff(1,1),fff(1,1),18*H_nElNode)
          enddo
        enddo  
      enddo
      call m_symLtoU (K(1,1),3*H_nElNode,3*H_nElNode)
      
      call Axv1_v2(ff(1,1),Beta(1),F(1),3*H_nElNode,3)
      call cv1_v1(F(1),3*H_nElNode,DH)
 
      return
      end