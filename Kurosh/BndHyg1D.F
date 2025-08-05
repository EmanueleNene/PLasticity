      subroutine BndHyg1D(BeamProps,Ang,Xi,Eps,MH)
c
c --- primary function: Finds longitudinal strains at end points of the beams in contact
c
c Input arguments:
c  BeamProps  (dp,arr(9,2),in)  -Beam properties as described below
c  Ang        (dp,sc,in)        -Contact angle
c  Xi         (dp,arr(2),in)    -Contact location (xi) on beam 1&2
c
c output arguments:
c  Eps      (dp,arr(2,2),out) -Strain at end poits of beam 1&2
c
      

c     & (L,beamprop(1)), (H,beamprop(2)), (W,beamprop(3)),
c     & (Ex,beamprop(4)), (Ez,beamprop(5)), (Gx,beamprop(6)),
c     & (alpx,beamprop(7)), (alpz,beamprop(8)), 
c     & (dh,beamprop(9))
      
c  Beam 1:
c   W1      (dp,sc,in)  -Width of cross-section 
c   H1      (dp,sc,in)  -Height of cross-section 
c   L1      (dp,sc,in)  -Beam Length
c   EL1     (dp,sc,in)  -Longitudianl elasticity modulus
c   ET1     (dp,sc,in)  -Transverse elasticity modulus
c   GLT1    (dp,sc,in)  -In-plane shear modulus (GLT)
c   NLT1    (dp,sc,in)  -In-plane Poisson's ratio
c   NTL1    (dp,sc,in)  -In-plane Poisson's ratio
c   BetaL1  (dp,sc,in)  -Longitudianl hydroexpansion coefficient
c   BetaT1  (dp,sc,in)  -Transverse hydroexpansion coefficient
c  Beam 2:
c   W2      (dp,sc,in)  -Width of cross-section 
c   H2      (dp,sc,in)  -Height of cross-section 
c   L2      (dp,sc,in)  -Beam Length
c   EL2     (dp,sc,in)  -Longitudianl elasticity modulus
c   ET2     (dp,sc,in)  -Transverse elasticity modulus
c   GLT2    (dp,sc,in)  -In-plane shear modulus (GLT)
c   NLT2    (dp,sc,in)  -In-plane Poisson's ratio
c   NTL2    (dp,sc,in)  -In-plane Poisson's ratio
c   BetaL2  (dp,sc,in)  -Longitudianl hydroexpansion coefficient
c   BetaT2  (dp,sc,in)  -Transverse hydroexpansion coefficient
c  Moisture
c   DH1     (dp,sc,in)  -Change in moisture content of beam 1
c   DH2     (dp,sc,in)  -Change in moisture content of beam 2
c

#include "impcom.inc"
      double precision    BeamProps(11,2),Ang,Xi(2)
      double precision    Eps(2,2), MH
      
      double precision   W1,H1,L1,EL1,ET1,GLT1,NLT1,NTL1,
     &                   BetaL1,BetaT1,DH1
      double precision   W2,H2,L2,EL2,ET2,GLT2,NLT2,NTL2,
     &                   BetaL2,BetaT2,DH2
      double precision   Ex,Betax,W01,W02,W0,W,D
      integer            i

      MH = 0.0d0
      call v_zero(Eps(1,1),4)
      
      call AssignProps(BeamProps,
     &                W1,H1,L1,EL1,ET1,GLT1,NLT1,NTL1,BetaL1,BetaT1,DH1,
     &                W2,H2,L2,EL2,ET2,GLT2,NLT2,NTL2,BetaL2,BetaT2,DH2)
      
      Ex = (cos(Ang))**4.0d0 * EL2 + (sin(Ang)**4.0d0) * ET2 + 
     &     4.0d0 * (sin(Ang))**2.0d0 * (cos(Ang))**2.0d0 *GLT2
      Betax = (cos(Ang))**2.0d0 * BetaL2 + (sin(Ang)**2.0d0) * BetaT2
      
      W0 = W2/sin(Ang)
      W01 = W0 * (1+BetaL1*DH1)
      W02 = W0 * (1+Betax *DH2)
      
      W = (EL1*H1 + Ex*H2) / (EL1*H1/W01 + Ex*H2/W02)
      
      D = W - W01
      
      Eps(1,1) = D/L1*( 4.0d0-3.0d0*(Xi(1)+1))
      Eps(2,1) = D/L1*(-2.0d0+3.0d0*(Xi(1)+1))
      
      
      
      Ex = (cos(Ang))**4.0d0 * EL1 + (sin(Ang)**4.0d0) * ET1 + 
     &     4.0d0 * (sin(Ang))**2.0d0 * (cos(Ang))**2.0d0 *GLT1
      Betax = (cos(Ang))**2.0d0 * BetaL1 + (sin(Ang)**2.0d0) * BetaT1
      
      W0 = W1/sin(Ang)
      W02 = W0 * (1+BetaL2*DH2)
      W01 = W0 * (1+Betax *DH1)
      
      W = (EL2*H2 + Ex*H1) / (EL2*H2/W02 + Ex*H1/W01)
      
      D = W - W02
      
      Eps(1,2) = D/L2*( 4.0d0-3.0d0*(Xi(2)+1))
      Eps(2,2) = D/L2*(-2.0d0+3.0d0*(Xi(2)+1))
      
      
      end subroutine BndHyg1D               