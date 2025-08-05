      subroutine AssignProps(BeamProps,
     &                W1,H1,L1,EL1,ET1,GLT1,NLT1,NTL1,BetaL1,BetaT1,DH1,
     &                W2,H2,L2,EL2,ET2,GLT2,NLT2,NTL2,BetaL2,BetaT2,DH2)
c --- primary function: Assign appropriate properties to variables
c
c Input arguments:
c  BeamProps  (dp,arr(9,2),in)  -Beam properties as described below
c
c output arguments:
c     & (L,beamprop(1)), (H,beamprop(2)), (W,beamprop(3)),
c     & (Ex,beamprop(4)), (Ez,beamprop(5)), (Gx,beamprop(6)),
c     & (alpx,beamprop(7)), (alpz,beamprop(8)), 
c     & (dh,beamprop(9))
c

      
#include "impcom.inc"

      double precision    BeamProps(11,2)
      
      double precision   W1,H1,L1,EL1,ET1,GLT1,NLT1,NTL1,
     &                   BetaL1,BetaT1,DH1
      double precision   W2,H2,L2,EL2,ET2,GLT2,NLT2,NTL2,
     &                   BetaL2,BetaT2,DH2
      
      double precision   ZERO
      parameter          (ZERO = 0.0d0) 
      
      L1      = BeamProps(1,1)
      H1      = BeamProps(2,1)
      W1      = BeamProps(3,1)
      EL1     = BeamProps(4,1)
      ET1     = BeamProps(5,1)
      GLT1    = BeamProps(6,1)
      BetaL1  = BeamProps(7,1)
      BetaT1  = BeamProps(8,1)
      DH1     = BeamProps(9,1)
      NLT1    = ZERO
      NTL1    = ZERO

      L2      = BeamProps(1,2)
      H2      = BeamProps(2,2)
      W2      = BeamProps(3,2)
      EL2     = BeamProps(4,2)
      ET2     = BeamProps(5,2)
      GLT2    = BeamProps(6,2)
      BetaL2  = BeamProps(7,2)
      BetaT2  = BeamProps(8,2)
      DH2     = BeamProps(9,2)
      NLT2    = ZERO
      NTL2    = ZERO
      
      return
      end