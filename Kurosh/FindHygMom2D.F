      subroutine FindHygMom2D(U,BeamProps,Ang,
     &                        KK,FF,Connectivity,ElemType,MH)
c --- primary function: Finds the in-plane moment caused by hygroexpansion

#include "impcom.inc"
#include "beam_hyg.inc"

      double precision    U(H_Dof),BeamProps(11,2),Ang
      double precision    KK(H_key(1)*H_nElNode,H_key(1)*H_nElNode,16)
      double precision    FF(H_key(1)*H_nElNode,16)
      integer             Connectivity(H_nElNode,H_nElem)
      integer             ElemType(H_nElem)

      double precision   MH
      
      
      double precision    U1(16),F1(16)
      integer             El, ElDof(16)
      
      integer             i,m,n,cc

      double precision   x0,y0,YX(16)
      double precision   dx1,dx2, dy
      double precision   W1,W2
      
      double precision   v_dot
      

      MH = 0.0d0
      W1 = BeamProps(3,1)
      W2 = BeamProps(3,2)

      dx1 =   W2/sin(Ang) / dble(H_key(3))
      dx2 =   W1/tan(Ang) / dble(H_key(3))
      dy  = W1 / dble(H_key(4))
      
      cc = H_key(3)*H_key(9)+H_key(7)
      
      y0 = -0.5d0* W1 - dy
      x0 = -0.5d0*( W2/sin(Ang)+ W1/tan(Ang) ) - dx2
      
      do m=1,H_key(4)
          y0 = y0 + dy
          x0 = x0 + dx2

          YX(1) = y0
          YX(3) = y0
          YX(5) = Y0 + dy
          YX(7) = y0 + dy
          if (H_Key(2) .eq. 2) then
              YX(9)  = y0
              YX(11) = y0 + 0.5d0*dy
              YX(13) = Y0 + dy
              YX(15) = y0 + 0.5d0*dy
          endif
          
          do n=1,H_key(3)
              El = cc + n
              YX(2) = x0 + (n-1)*dx1
              YX(4) = YX(2) + dx1
              YX(6) = YX(4) + dx2
              YX(8) = YX(2) + dx2
              if (H_Key(2) .eq. 2) then
                  YX(10) = YX(2) + 0.5d0*dx1
                  YX(12) = YX(4) + 0.5d0*dx2
                  YX(14) = YX(8) + 0.5d0*dx1
                  YX(16) = YX(2) + 0.5d0*dx2
              endif
              
              do i = 1, H_nElNode
                  ElDof(2*i-1) = Connectivity(i,El)*2-1
                  ElDof(2*i  ) = Connectivity(i,El)*2
              enddo
              
              U1(1:2*H_nElNode) = U(ElDof(1:2*H_nElNode))
              
              call Axv1_v2 (KK(1,1,ElemType(El)), U1(1), F1(1),
     &                      2*H_nElNode,2*H_nElNode)

c Hygro deformations are ignored as they are symmetric
c This part should be modified to include the hollow case!!!!
              MH = MH + v_dot(F1,YX,2*H_nElNode)

              cc = cc+1
          enddo
          cc = cc+H_key(7)+H_key(3)+H_key(8)
      enddo

      
      return
      end
          
      