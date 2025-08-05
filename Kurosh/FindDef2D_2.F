      subroutine FindDef2D_2(U,BeamProps,Ang,Connectivity,BndEps,W)
c --- primary function: Finds deformation at the middle of the bond 

#include "impcom.inc"
#include "beam_hyg.inc"

      double precision    U(H_Dof),BeamProps(11,2),Ang
      integer             Connectivity(H_nElNode,H_nElem)
      double precision    BndEps(8,2),W(2)
      
      double precision    J(2,2),Jinv(2,2),detJ
      double precision    Bxi(2,H_nElNode)
c      double precision    B(3,2*H_nElNode)
      double precision    B1(2*H_nElNode,H_nGpt*H_nGpt)
      double precision    B2(2*H_nElNode,H_nGpt*H_nGpt)
      
      double precision    XY(H_nElNode,2),XY2(H_nElNode,2)
      
      double precision    U1(2*H_nElNode),U2(2*H_nElNode)
      integer             El, ElDof(2*H_nElNode)
      
      double precision    EpsX1, EpsX2
      
      double precision    nis, sss
      
      double precision     Hlf
      parameter           (Hlf    = 0.5d0)
      
      double precision    v_dot
      
      integer             i,p,m,n,cc, count

      double precision   W1,W2, BetaL1,BetaL2,DH1,DH2
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
      
      count = 1
      do i=1,H_nGpt
          do p=1,H_nGpt
             call B_xi_2(H_xG(i),H_xG(p),H_nElNode,Bxi(1,1))
             
             call JacH(Bxi(1,1),XY(1,1),H_key(1),H_nElNode,
     &                J(1,1),Jinv(1,1),detJ)
             call B_Fcn_2x(Bxi(1,1),Jinv(1,1),H_nElNode,B1(1,count))
             
             call JacH(Bxi(1,1),XY2(1,1),H_key(1),H_nElNode,
     &                J(1,1),Jinv(1,1),detJ)
             call B_Fcn_2x(Bxi(1,1),Jinv(1,1),H_nElNode,B2(1,count))
             
             count = count + 1
          enddo
      enddo
      
      EpsX1 = 0.0d0
      EpsX2 = 0.0d0
      
      cc = H_key(3)*H_key(9)+H_key(7)
      do m=1,H_key(4)
          do n=1,H_key(3)
              El = cc + n
              
              do i = 1, H_nElNode
                  ElDof(2*i-1) = Connectivity(i,El)*2-1
                  ElDof(2*i  ) = Connectivity(i,El)*2
              enddo
              
              U1(1:2*H_nElNode) = U(ElDof(1:2*H_nElNode))
              do i = 1, H_nElNode
                  call Axv1_v2(Rot(1,1),U1(i*2-1),U2(i*2-1),2,2)
              enddo
              
              count = 1
              do i=1,H_nGpt
                  do p=1,H_nGpt
                      sss = H_wG(i)*H_wG(p)
c                      sss2 = detJ*wx(i)*wy(p)
                      EpsX1 = EpsX1 + 
     &                         sss  * v_dot(B1(1,count),U1,2*H_nElNode)
                      EpsX2 = EpsX2 +
     &                         sss  * v_dot(B2(1,count),U2,2*H_nElNode)
                      count = count +1
                  enddo
              enddo
              cc = cc+1
          enddo
          cc = cc+H_key(7)+H_key(3)+H_key(8)
      enddo

      nis  = 1.0d0 / sin(Ang)
      
      EpsX1 = 0.25d0 * EpsX1/(H_key(3)*H_key(4)) - BetaL1*DH1
      EpsX2 = 0.25d0 * EpsX2/(H_key(3)*H_key(4)) - BetaL2*DH2
      
      W(1) = (1 + EpsX2+BetaL2*DH2) * W1 * nis
      W(2) = (1 + EpsX1+BetaL1*DH1) * W2 * nis
      
      BndEps(:,1) = EpsX1
      BndEps(:,2) = EpsX2
      
      return
      end
          
      