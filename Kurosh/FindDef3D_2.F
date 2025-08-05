      subroutine FindDef3D_2(U,BeamProps,ac,bc,dir,ta,Ang,
     &                       Connectivity,BndEps)
c --- primary function: Finds deformation at the middle of the bond 

#include "impcom.inc"
#include "beam_hyg.inc"

      double precision    U(H_Dof),BeamProps(11,2),Ang
      double precision    ac(3),bc(3),dir(3),ta(3),tb(3)
      integer             Connectivity(H_nElNode,H_nElem)
      double precision    BndEps(8,2)
      
      double precision    J(3,3),Jinv(3,3),detJ
      double precision    Bxi(3,H_nElNode)
c      double precision    B(3,2*H_nElNode)
      double precision    B1(3*H_nElNode,8)
      double precision    B2(3*H_nElNode,8)
      double precision    BB(3*H_nElNode)
      
      double precision    XYc(4,2), XYc2(4,2)
      double precision    XY(H_nElNode,3),XY2(H_nElNode,3)
      
      double precision    U1(3*H_nElNode),U2(3*H_nElNode),V(3*H_nElNode)
      double precision    UU(3*H_nElNode)
      integer             El, ElDof(3*H_nElNode)
      
      double precision    GpBndEps(8,2)
      
      double precision    nis, sss
      
      double precision     Hlf
      parameter           (Hlf    = 0.5d0)
      
      double precision    v_dot
      
      integer             i,iy,m,BM, count

      double precision   H1,H2,W1,W2
      double precision   BetaL(2),DH
      double precision   Rot(3,3),RotT(2,2)
      double precision   c,s
      
      double precision   Gx,x,y,z,ss
      integer            nn(5,2),PtOrder(3,8)
      integer            nxyz(3,2),nElem(2,2,2),nEl0
      double precision   sxyz(3,2)
      
      double precision    dXYZ_L(3),XYZ_G(3), dXYZ_G(3,3)
      double precision    XYZ(3),mst
      
      
      call NdIOrder(PtOrder)
      Gx     = 0.5773502691896258

      H1      = BeamProps(2,1)/H_key(11)
      W1      = BeamProps(3,1)/H_key(4)
      BetaL(1)= BeamProps(7,1)
      H2      = BeamProps(2,2)/H_key(12)
      W2      = BeamProps(3,2)/H_key(3)
      BetaL(2)= BeamProps(7,2)
      c = cos(Ang)
      s = sin(Ang)
      call v_zero (Rot(1,1),9)
      Rot(1,1) =  c
      Rot(1,2) =  s
      Rot(2,1) = -s
      Rot(2,2) =  c
      Rot(3,3) =  1
      RotT(1,1) =  c
      RotT(1,2) = -s
      RotT(2,1) =  s
      RotT(2,2) =  c
      
      call BndPts(W1,W2,Ang,XYc(1:4,:))
      call AxB_C(XYc(1:4,:),RotT(1,1),XYc2(1:4,:),4,2,4,4,2,2)

      XY(1:4,1:2) = XYc
      XY(1:4,3)   = 0.0d0
      XY(5:8,1:2) = XYc
      XY(5:8,3)   =  H1
      if (H_nElNode .eq. 20) then
          XY(9 ,1:2) = Hlf * (XYc(2,:)+XYc(3,:))
          XY(10,1:2) = Hlf * (XYc(3,:)+XYc(4,:))
          XY(11,1:2) = Hlf * (XYc(4,:)+XYc(1,:))
          XY(12,1:2) = Hlf * (XYc(1,:)+XYc(2,:))
          XY(9:12,3) = 0.0d0
          XY(13:16,1:2) = XY(9:12,1:2)
          XY(13:16,3)   =  H1
          XY(17:19,1:2) = XYc(2:4,:)
          XY(20,1:2)    = XYc(1,:)
          XY(17:20,3)   = Hlf * H1
      endif
      
      XY2(1:4,1:2) = XYc2
      XY2(1:4,3)   = 0.0d0
      XY2(5:8,1:2) = XYc2
      XY2(5:8,3)   =  H2
      if (H_nElNode .eq. 20) then
          XY2(9 ,1:2) = Hlf * (XYc2(2,:)+XYc2(3,:))
          XY2(10,1:2) = Hlf * (XYc2(3,:)+XYc2(4,:))
          XY2(11,1:2) = Hlf * (XYc2(4,:)+XYc2(1,:))
          XY2(12,1:2) = Hlf * (XYc2(1,:)+XYc2(2,:))
          XY2(9:12,3) = 0.0d0
          XY2(13:16,1:2) = XY2(9:12,1:2)
          XY2(13:16,3)   =  H2
          XY2(17:19,1:2) = XYc2(2:4,:)
          XY2(20,1:2)    = XYc2(1,:)
          XY2(17:20,3)   = Hlf * H2
      endif
      
      nn(1,1) = H_key(3)
      nn(2,1) = H_key(4)
      nn(3,1) = H_key(11)
      nn(4,1) = H_key(9)
      nn(5,1) = H_key(10)
      
      nn(1,2) = H_key(4)
      nn(2,2) = H_key(3)
      nn(3,2) = H_key(12)
      nn(4,2) = H_key(7)
      nn(5,2) = H_key(8)
      
      nEl0 = 0
      do BM=1,2
          do i=1,3
              ss = (1.0d0-Gx)*nn(i,BM)/2.0d0
              nxyz(i,1) = int(ss) + 1
              sxyz(i,1) = 2.0d0 * (ss - dble(nxyz(i,1)) ) + 1.0d0
              
              ss = (1.0d0+Gx)*nn(i,BM)/2.0d0
              nxyz(i,2) = int(ss) + 1
              sxyz(i,2) = 2.0d0 * (ss - dble(nxyz(i,2)) ) + 1.0d0
          enddo
          
          iy = 3-BM
          nElem(1,BM,1)   = nEl0 +
     &               nn(1,BM)*(nn(2,BM)+nn(4,BM)+nn(5,BM))*(nxyz(3,1)-1)
     &              +(nn(4,BM)+nxyz(2,1)-1)*nn(1,BM) + nxyz(1,1)
          nElem(2,BM,1) = nElem(1,BM,1) + nxyz(1,2)-nxyz(1,1)
          nElem(:,iy,1) = nElem(:,BM,1) + (nxyz(2,2)-nxyz(2,1))*nn(1,BM)
          
          nElem(:,:,2) = nElem(:,:,1) +
     &       nn(1,BM)*(nn(2,BM)+nn(4,BM)+nn(5,BM))*(nxyz(3,2)-nxyz(3,1))
          
          nEl0 = nn(1,BM)*(nn(2,BM)+nn(4,BM)+nn(5,BM))*nn(3,BM)
          
          do count=1,8
              x = sxyz(1,PtOrder(1,count))
              y = sxyz(2,PtOrder(2,count))
              z = sxyz(3,PtOrder(3,count))
              
            El=nElem(PtOrder(1,count),PtOrder(2,count),PtOrder(3,count))
              
              

              call B_xi_3(x,y,z,H_nElNode,Bxi(1,1))
              call JacH(Bxi(1,1),XY(1,1),H_key(1),H_nElNode,
     &                    J(1,1),Jinv(1,1),detJ)
              call B_Fcn_3x(Bxi(1,1),Jinv(1,1),H_nElNode,BB(1))

              do m = 1, H_nElNode
                  ElDof(3*m-2) = Connectivity(m,El)*3-2
                  ElDof(3*m-1) = Connectivity(m,El)*3-1
                  ElDof(3*m  ) = Connectivity(m,El)*3
              enddo
              
              UU(1:3*H_nElNode) = U(ElDof(1:3*H_nElNode))
              
              if (BM .eq. 2) then
                  do i = 1, H_nElNode
                      call Axv1_v1(Rot(1,1),UU(i*3-2),3,3)
                  enddo
              endif
              
              
              if ( ( H_key(16) .eq. 0) .or.
     &             (BeamProps(11,BM) .eq. 0.0d0) ) then
                  DH = BeamProps(9,BM)
              else

                  call ElemLocX(El,BeamProps,Ang,dXYZ_L)
                  call LocToGlobal(El,dXYZ_L,ac,bc,dir,ta,
     &                             XYZ_G,dXYZ_G)
                  
                  XYZ = XYZ_G + x*dXYZ_G(:,1)
     &                        + y*dXYZ_G(:,2)
     &                        + z*dXYZ_G(:,3)
                  
                  call GetTabMoist(BeamProps(11,BM),XYZ(1),mst)
                  DH = mst - BeamProps(10,BM)

c                  if (abs(DH-BeamProps(9,BM)).gt.1.4d0) then
c                      write (6,*)'DH',DH,BeamProps(9,BM)
c                      write (6,*)'Height',DH,BeamProps(2,BM)
c                      write (6,*)'dir',dir
c                      write (6,*)'ta',ta
c                      write (6,*)'x,y,z',x,y,z
c                      write (6,*)'XYZ',XYZ
c                      write (6,*)'XYZ_G',XYZ_G
c                      write (6,*)'dXYZ_L',dXYZ_L
c                      write (6,*)'dXYZ_G',dXYZ_G
c                  endif
                      
              endif
              
              
              GpBndEps(count,BM) =v_dot(BB(1),UU,3*H_nElNode)
     &                                - BetaL(BM)*DH
              
          enddo
          XY = XY2
  
      enddo

      call GpBndEpsToNd(GpBndEps,BndEps)
      
      return
      end
          
      