      subroutine SetupKKandFF(BeamProps,Ang,KK,FF)
      
#include "impcom.inc"
#include "beam_hyg.inc"
      
      double precision    BeamProps(11,2),Ang
     
      double precision    KK(H_key(1)*H_nElNode,H_key(1)*H_nElNode,
     &                                                     H_ElemTypK)
      double precision    FF(H_key(1)*H_nElNode,H_ElemTypF)
      
      double precision   W1,H1,L1,EL1,ET1,GLT1,NLT1,NTL1,
     &                   BetaL1,BetaT1,DH1
      double precision   W2,H2,L2,EL2,ET2,GLT2,NLT2,NTL2,
     &                   BetaL2,BetaT2,DH2
      
      double precision    X(4,2), X2(4,2)
      
      double precision    Rot(3,3),RotT(3,3)
      double precision    T(60,60), TT(60,60)
      
      
      call AssignProps(BeamProps,
     &                W1,H1,L1,EL1,ET1,GLT1,NLT1,NTL1,BetaL1,BetaT1,DH1,
     &                W2,H2,L2,EL2,ET2,GLT2,NLT2,NTL2,BetaL2,BetaT2,DH2)
      
      call RotMtx(Ang,Rot,RotT,T,TT,H_key(1),H_key(1)*H_nElNode)
      
      if (H_key(1) .eq. 2) then
          if ((H_key(5) .eq. 0) .and. (H_key(6) .eq. 0)) then               both beams have solid sections
              
              call BndPts(W1/H_key(4),W2/H_key(3),Ang,X(1,1))
              
              call Calc_K_F_2(X,EL1,ET1,GLT1,NLT1,NTL1,H1,BetaL1,BetaT1,
     &                      DH1,KK(:,:,1),FF(:,1))
              
              call AxB_C(X(1,1),RotT(1,1),X2(1,1),4,3,4,4,2,2)
              call Calc_K_F_2(X2,EL2,ET2,GLT2,NLT2,NTL2,H2,BetaL2,BetaT2
     &                     ,DH2,KK(:,:,2),FF(:,2))
              call BtxAxB_A(KK(:,:,2),2*H_nElNode,T(1,1),60,
     &                      2*H_nElNode,2*H_nElNode)
              call Axv1_v1(TT(1,1),FF(:,2),60,2*H_nElNode)
              
              KK(:,:,9)  = KK(:,:,1)
              KK(:,:,11) = KK(:,:,1)
              KK(:,:,13) = KK(:,:,2)
              KK(:,:,15) = KK(:,:,2)
              
              FF(:,9)  = FF(:,1)
              FF(:,11) = FF(:,1)
              FF(:,13) = FF(:,2)
              FF(:,15) = FF(:,2)

          elseif (H_key(5) .eq. 0) then
              
          elseif (H_key(6) .eq. 0) then
              
          else
              
          endif
      elseif (H_key(1) .eq. 3) then
          if (H_key(16) .eq. 2) then
              return
          endif
          
          if ((H_key(5) .eq. 0) .and. (H_key(6) .eq. 0)) then               both beams have solid sections
              call BndPts(W1/H_key(4),W2/H_key(3),Ang,X(1,1))
              
              call Calc_K_F_3(X,EL1,ET1,GLT1,NLT1,NTL1,H1/H_key(11),
     &                        BetaL1,BetaT1,DH1,KK(:,:,1),FF(:,1))
              
              call AxB_C(X(1,1),RotT(1,1),X2(1,1),4,3,4,4,2,2)
              call Calc_K_F_3(X2,EL2,ET2,GLT2,NLT2,NTL2,H2/H_key(12),
     &                        BetaL2,BetaT2,DH2,KK(:,:,2),FF(:,2))
              call BtxAxB_A(KK(:,:,2),3*H_nElNode,T(1,1),60,
     &                      3*H_nElNode,3*H_nElNode)
              call Axv1_v1(TT(1,1),FF(:,2),60,3*H_nElNode)
              
          endif
      endif
      
      return
      end