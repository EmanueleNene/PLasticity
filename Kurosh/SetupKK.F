      subroutine SetupKK(Ang,BB,WW,CC,KK)
      
#include "impcom.inc"
#include "beam_hyg.inc"
      
      double precision    BB(6,3*H_nElNode,H_nGpt*H_nGpt*H_nGpt,
     &                                                     H_ElemTypB)
      double precision    WW(H_nGpt*H_nGpt*H_nGpt,H_ElemTypB)
      double precision    CC(6,6,2)
     
      double precision    KK(H_key(1)*H_nElNode,H_key(1)*H_nElNode,
     &                                                     H_ElemTypK)
      
      double precision    fff(3*H_nElNode,6)
     
      double precision    Ang
      double precision    Rot(3,3),RotT(3,3)
      double precision    T(60,60), TT(60,60)
      
      integer             i,p,q,count
      
      
      call RotMtx(Ang,Rot,RotT,T,TT,H_key(1),H_key(1)*H_nElNode)
      
      if (H_key(1) .eq. 3) then
          
          if ((H_key(5) .eq. 0) .and. (H_key(6) .eq. 0)) then               both beams have solid sections
              call v_zero (KK(1,1,1),9*H_nElNode*H_nElNode*H_ElemTypK)
      
              count = 0
              do i=1,H_nGpt
                do p=1,H_nGpt
                  do q=1,H_nGpt
                      count = count+1

                     call cAtxBxA_C(BB(1,1,count,1),CC(1,1,1),KK(1,1,1),
     &               6,6,3*H_nElNode,6,3*H_nElNode,fff(1,1),WW(count,1))
                     
                     call cAtxBxA_C(BB(1,1,count,2),CC(1,1,2),KK(1,1,2),
     &               6,6,3*H_nElNode,6,3*H_nElNode,fff(1,1),WW(count,2))
                     
                  enddo
                enddo  
              enddo
              call m_symLtoU (KK(1,1,1),3*H_nElNode,3*H_nElNode)
              call m_symLtoU (KK(1,1,2),3*H_nElNode,3*H_nElNode)
              
              call BtxAxB_A(KK(1,1,2),3*H_nElNode,T(1,1),60,
     &                      3*H_nElNode,3*H_nElNode)
              
          endif
      endif
      
      return
      end