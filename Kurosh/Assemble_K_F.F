      subroutine Assemble_K_F(BeamProps,ac,bc,dir,ta,Ang,KK,FF,BB,WW,
     &                        CC,Connectivity,ElemType,K,F)
      
#include "impcom.inc"
#include "beam_hyg.inc"

      double precision    BeamProps(11,2),Ang
      double precision    ac(3),bc(3),dir(3),ta(3),tb(3)

      double precision    K(H_Dof,H_Dof)
      double precision    F(H_Dof)
      double precision    KK(H_key(1)*H_nElNode,H_key(1)*H_nElNode,
     &                                                     H_ElemTypK)
      double precision    FF(H_key(1)*H_nElNode,H_ElemTypF)
      double precision    BB(6,3*H_nElNode,H_nGpt*H_nGpt*H_nGpt,
     &                                                     H_ElemTypB)
      double precision    WW(H_nGpt*H_nGpt*H_nGpt,H_ElemTypB)
      double precision    CC(6,6,2)
      
      integer             Connectivity(H_nElNode,H_nElem)
      integer             ElemType(H_nElem)
      
      integer             i, j
      integer             loci(H_key(1)*H_nElNode)
      
      double precision    FFF(H_key(1)*H_nElNode)
      
      call v_zero(K(1,1),H_Dof*H_Dof)
      call v_zero(F(1),H_Dof)
      
      do i=1,H_nElem
          if (H_key(1) .eq. 2) then
              do j=1,H_nElNode
                  loci(j*2-1) = Connectivity(j,i)*2-1
                  loci(j*2)   = Connectivity(j,i)*2
              enddo
          elseif (H_key(1) .eq. 3) then
              do j=1,H_nElNode
                  loci(j*3-2) = Connectivity(j,i)*3-2
                  loci(j*3-1) = Connectivity(j,i)*3-1
                  loci(j*3)   = Connectivity(j,i)*3
              enddo
          endif
          
          K(loci,loci) = K(loci,loci) + KK(:,:,ElemType(i))
          if (H_key(16) .eq. 0) then
              F(loci)      = F(loci)      + FF(:  ,ElemType(i))
          else
              call Calc_F_3(BeamProps,ac,bc,dir,ta,ElemType(i),
     &                                       i,Ang,BB,WW,CC,FFF)
              F(loci)      = F(loci)      + FFF
          endif
              
      enddo
      
      return
      end
