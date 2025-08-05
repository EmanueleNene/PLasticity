      subroutine Calc_F_3(BeamProps,ac,bc,dir,ta,Mat,ElNum,Ang,BB,WW,
     &                    CC,F)
      
#include "impcom.inc"
#include "beam_hyg.inc"
      double precision    BeamProps(11,2),Ang
      double precision    ac(3),bc(3),dir(3),ta(3),tb(3)
      integer             Mat, ElNum
      double precision    BB(6,3*H_nElNode,H_nGpt*H_nGpt*H_nGpt,
     &                                                     H_ElemTypB)
      double precision    WW(H_nGpt*H_nGpt*H_nGpt,H_ElemTypB)
      double precision    CC(6,6,2)
     
      double precision    F(3*H_nElNode)
      
      double precision    ff(3*H_nElNode,6),fff(3*H_nElNode,6)
      
      double precision    Beta(6), DH, mst
      double precision    dXYZ_L(3),XYZ_G(3), dXYZ_G(3,3)
      double precision    XYZ(3)
      
      double precision    Rot(3,3),RotT(3,3)
      double precision    T(60,60), TT(60,60)
      
      integer             i,p,q,count
      
      
      if (H_key(1) .eq. 3) then
          
          Beta(1) = BeamProps(7,Mat)
          Beta(2) = BeamProps(8,Mat)
          Beta(3) = BeamProps(8,Mat)
          Beta(4) = 0.0d0
          Beta(5) = 0.0d0
          Beta(6) = 0.0d0
          
          if ((H_key(5) .eq. 0) .and. (H_key(6) .eq. 0)) then               both beams have solid sections
              call v_zero (F(1),3*H_nElNode)
              call v_zero (ff(1,1),18*H_nElNode)
              call v_zero (fff(1,1),18*H_nElNode)
              
              call ElemLocX(ElNum,BeamProps,Ang,dXYZ_L)
              call LocToGlobal(ElNum,dXYZ_L,ac,bc,dir,ta,
     &                         XYZ_G,dXYZ_G)
              count = 0

              do i=1,H_nGpt
                do p=1,H_nGpt
                  do q=1,H_nGpt
                      count = count+1
                      XYZ = XYZ_G + H_xG(i)*dXYZ_G(:,1)
     &                            + H_xG(p)*dXYZ_G(:,2)
     &                            + H_xG(q)*dXYZ_G(:,3)
                      
                      if (BeamProps(11,Mat) .eq. 0.0d0) then
                          DH = BeamProps(9,Mat)
                      else
                          call GetTabMoist(BeamProps(11,Mat),XYZ(1),mst)
                          DH = mst - BeamProps(10,Mat)
                      endif

                      fff = WW(count,Mat)*DH*
     &                  matmul(transpose(BB(:,:,count,Mat)),CC(:,:,Mat))
                      
                     call v1pv2_v1(ff(1,1),fff(1,1),18*H_nElNode)
                     
                  enddo
                enddo  
              enddo
              
              call Axv1_v2(ff(1,1),Beta(1),F(1),3*H_nElNode,3)
              
              if (Mat .eq. 2) then
                  call RotMtx(Ang,Rot,RotT,T,TT,3,3*H_nElNode)
                  F = Matmul(TT(1:3*H_nElNode,1:3*H_nElNode),F)
c                  call Axv1_v1(TT(1,1),F(1),3*H_nElNode,3*H_nElNode)
              endif
              
          endif
      endif
      
      return
      end