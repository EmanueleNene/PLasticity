      subroutine FindDef2D_1(U,BeamProps,Ang,Connectivity,BndEps)
c --- primary function: Finds deformation at the middle of the bond 

#include "impcom.inc"
#include "beam_hyg.inc"

      double precision    U(H_Dof),BeamProps(11,2),Ang
      integer             Connectivity(H_nElNode,H_nElem)
      double precision    BndEps(8,2)
      
      double precision    J(2,2),Jinv(2,2),detJ
      double precision    Bxi(2,H_nElNode)
c      double precision    B(3,2*H_nElNode)
      double precision    B1(2*H_nElNode),B2(2*H_nElNode)
      
      double precision    XY(H_nElNode,2),XY2(H_nElNode,2)
      
      double precision    U1(16),U2(16)
      integer             El, ElDof(16)
      
      double precision    x,y
      
      double precision    EpsX1, EpsX2
      
c      double precision    nis
      
      double precision     Hlf
      parameter           (Hlf    = 0.5d0)
      
      double precision    v_dot
      
      integer             i
      
      double precision   W1,W2,BetaL1,BetaL2,DH1,DH2
      double precision   Rot(2,2),RotT(2,2)
      double precision   c,s

      W1      = BeamProps(3,1)
      BetaL1  = BeamProps(7,1)
      DH1     = BeamProps(9,1)
      W2      = BeamProps(3,2)
      BetaL2  = BeamProps(7,2)
      DH2     = BeamProps(9,2)
      c = cos(Ang)
      s = sin(Ang)
      Rot(1,1) =  c
      Rot(1,2) =  s
      Rot(2,1) = -s
      Rot(2,2) =  c
      RotT(1,1) =  c
      RotT(1,2) = -s
      RotT(2,1) =  s
      RotT(2,2) =  c
      
      call BndPts(W1/H_key(4),W2/H_key(3),Ang,XY(1:4,:))
      call AxB_C(XY(1:4,:),RotT(1,1),XY2(1:4,:),4,2,4,4,2,2)

      if (H_nElNode .eq. 8) then
          XY(5,:) = Hlf * (XY(1,:)+XY(2,:))
          XY(6,:) = Hlf * (XY(2,:)+XY(3,:))
          XY(7,:) = Hlf * (XY(3,:)+XY(4,:))
          XY(8,:) = Hlf * (XY(4,:)+XY(1,:))
          
          XY2(5,:) = Hlf * (XY2(1,:)+XY2(2,:))
          XY2(6,:) = Hlf * (XY2(2,:)+XY2(3,:))
          XY2(7,:) = Hlf * (XY2(3,:)+XY2(4,:))
          XY2(8,:) = Hlf * (XY2(4,:)+XY2(1,:))
      endif
      
      El = H_key(3)*H_key(9) + 
     &    (H_key(4)-1)/2*(2*H_key(3)+H_key(7)+H_key(8)) +
     &     H_key(7)+H_key(3)
      
      if (   ( (H_key(3)/2)*2 .eq. H_key(3) ) .and.
     &       ( (H_key(4)/2)*2 .eq. H_key(4) )   ) then
          x = 1.0d0
          y = 1.0d0
      elseif ( (H_key(3)/2)*2 .eq. H_key(3) ) then
          x = 1.0d0
          y = 0.0d0
      elseif ( (H_key(4)/2)*2 .eq. H_key(4) ) then
          x = 0.0d0
          y = 1.0d0
      else
          x = 0.0d0
          y = 0.0d0
      endif
      
      call B_xi_2(x,y,H_nElNode,Bxi(1,1))
      
      call JacH(Bxi(1,1), XY(1,1),H_key(1),H_nElNode,
     &         J(1,1),Jinv(1,1),detJ)
      call B_Fcn_2x(Bxi(1,1),Jinv(1,1),H_nElNode,B1(1))
      
      call JacH(Bxi(1,1), XY2(1,1),H_key(1),H_nElNode,
     &         J(1,1),Jinv(1,1),detJ)
      call B_Fcn_2x(Bxi(1,1),Jinv(1,1),H_nElNode,B2(1))
      
      do i = 1, H_nElNode
          ElDof(2*i-1) = Connectivity(i,El)*2-1
          ElDof(2*i  ) = Connectivity(i,El)*2
      enddo
      
      U1(1:2*H_nElNode) = U(ElDof(1:2*H_nElNode))
      do i = 1, H_nElNode
          call Axv1_v2(Rot(1,1),U1(i*2-1),U2(i*2-1),2,2)
      enddo

      EpsX1 = v_dot(B1(1),U1,2*H_nElNode) - BetaL1*DH1
      EpsX2 = v_dot(B2(1),U2,2*H_nElNode) - BetaL2*DH2
      
      
      BndEps(:,1) = EpsX1
      BndEps(:,2) = EpsX2
c      nis  = 1.0d0 / sin(Ang)
c      Def(1) = EpsX1 * W2 * nis
c      Def(2) = EpsX2 * W1 * nis
      
      return
      end
          
      