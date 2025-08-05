      subroutine C_Mtx(EL,ET,GLT,GTT,NLT,NTL,NTT,CDim,C)
c --- primary function: Calculate material stiffness matrix
c
c Input arguments:
c   EL     (dp,sc,in)  -Longitudianl elasticity modulus
c   ET     (dp,sc,in)  -Transverse elasticity modulus
c   GLT    (dp,sc,in)  -In-plane shear modulus (GLT)
c   GTT    (dp,sc,in)  -Cross-section shear modulus (GTT)
c   NLT    (dp,sc,in)  -In-plane Poisson's ratio
c   NTL    (dp,sc,in)  -In-plane Poisson's ratio
c   NTT    (dp,sc,in)  -Cross-section Poisson's ratio
c  CDim    (int,sc,in) - Dimension of C-matrix (3 for 2D and 6 for 3D)
c
c output arguments:
c  C     (dp,arr(CDim,CDim),out) - Material stiffness matrix
c
#include "impcom.inc"
c
 
      
      integer             CDim
      double precision    EL,ET,GLT,GTT,NLT,NTL,NTT
      double precision    C(CDim,CDim)
      
      double precision    s,sLT,sTT
      
      call v_zero (C(1,1),CDim*CDim)

      
      if (CDim .eq. 3) then                                             (Plane stress)
          s = 1.0d0 / (1.0d0 - NLT*NTL)
          
          C(1,1) = EL * s
          C(1,2) = EL * s * NTL
          C(2,1) = ET * s * NLT
          C(2,2) = ET * s
          C(3,3) = GLT
      elseif (CDim .eq. 6) then
          sLT = NLT*NTL
          sTT = NTT*NTT
          s   = 1.0d0 / (1.0d0 -2.0d0*sLT -sTT -2.0d0*sLT*NTT)
          
          C(1,1) = EL * s * (1.0d0-sTT)
          C(1,2) = EL * s * (NTL+NTL*NTT)
          C(1,3) = C(1,2)
          C(2,1) = ET * s * (NLT+NLT*NTT)
          C(2,2) = ET * s * (1.0d0-sLT)
          C(2,3) = ET * s * (NTT+sLT)
          C(3,1) = C(2,1)
          C(3,2) = C(2,3)
          C(3,3) = C(2,2)
          C(4,4) = 2.0d0*GLT
          C(5,5) = 2.0d0*GLT
          C(6,6) = 2.0d0*GTT
      else
          call erhandler ('BndHygStrn>C_Mtx',5000,2,
     &        'Wrong C-matrix dimension',0.0d0,' ')
      endif
     
999   return
      end      