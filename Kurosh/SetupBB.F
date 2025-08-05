      subroutine SetupBB(BeamProps,Ang,BB,WW,CC)
      
#include "impcom.inc"
#include "beam_hyg.inc"
      
      double precision    BeamProps(11,2),Ang
     
      double precision    BB(6,3*H_nElNode,H_nGpt*H_nGpt*H_nGpt,
     &                                                     H_ElemTypB)
      double precision    WW(H_nGpt*H_nGpt*H_nGpt,H_ElemTypB)
      double precision    CC(6,6,2)
      
      double precision   W1,H1,L1,EL1,ET1,GLT1,NLT1,NTL1,
     &                   BetaL1,BetaT1,DH1
      double precision   W2,H2,L2,EL2,ET2,GLT2,NLT2,NTL2,
     &                   BetaL2,BetaT2,DH2
      
      double precision   GTT,NTT
      
      double precision    X(4,2), X2(4,2), c,s,RotT(2,2)
      
      
      call AssignProps(BeamProps,
     &                W1,H1,L1,EL1,ET1,GLT1,NLT1,NTL1,BetaL1,BetaT1,DH1,
     &                W2,H2,L2,EL2,ET2,GLT2,NLT2,NTL2,BetaL2,BetaT2,DH2)
      
      c = cos(Ang)
      s = sin(Ang)
      RotT(1,1) =  c
      RotT(1,2) = -s
      RotT(2,1) =  s
      RotT(2,2) =  c
      
      if (H_key(1) .eq. 3) then
          if ((H_key(5) .eq. 0) .and. (H_key(6) .eq. 0)) then               both beams have solid sections
              
              NTT = 0.0d0
              GTT = GLT1*0.1d0
              call C_Mtx(EL1,ET1,GLT1,GTT,NLT1,NTL1,NTT,6,CC(1,1,1))
              GTT = GLT2*0.1d0
              call C_Mtx(EL2,ET2,GLT2,GTT,NLT2,NTL2,NTT,6,CC(1,1,2))
              
              
              call BndPts(W1/H_key(4),W2/H_key(3),Ang,X(1,1))
              
              call Calc_B_3(X,H1/H_key(11),BB(1,1,1,1),WW(1,1))
              

              call AxB_C(X(1,1),RotT(1,1),X2(1,1),4,2,4,4,2,2)
              
              call Calc_B_3(X2,H2/H_key(12),BB(1,1,1,2),WW(1,2))
              
          endif
      endif
      
      return
      end