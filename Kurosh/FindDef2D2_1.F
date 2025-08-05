      subroutine FindDef2D2_1(U,Ang,Def)
c --- primary function: Finds deformation at the middle of the bond 
     
#include "impcom.inc"
#include "beam_hyg.inc"

      double precision    U(H_Dof),Ang
      double precision    Def(2)
      
      integer             N1, N2, N3, N4
      double precision    U1, U2
      double precision    c , s

      c = cos(Ang)
      s = sin(Ang)
      
      if (2*(H_key(4)/2) .eq. H_key(4)) then
          N1 = (3*H_key(3)+2) * H_key(9) +
     &         (H_key(4)/2) *( 3*(H_key(7)+H_key(3)+H_key(8))+2 ) +
     &         2*H_key(7) +1
          N2 = N1 + 2*H_key(3)
          Def(1) = U(N2*2-1) - U(N1*2-1)
      else
          N1 = (3*H_key(3) +2) * H_key(9) +
     &         (H_key(4)/2 +1) *( 3*(H_key(7)+H_key(3)+H_key(8))+2 ) -
     &         (H_key(3)+H_key(8))
          N2 = N1 + H_key(3)
          Def(1) = U(N2*2-1) - U(N1*2-1)
      endif
      
      
          N1 = (3*H_key(3)+2) * H_key(9) +
     &         2*H_key(7) +H_key(3) + 1
          N2 = (3*H_key(3)+2) * H_key(9) +
     &         H_key(4) *( 3*(H_key(7)+H_key(3)+H_key(8))+2 ) +
     &         2*H_key(7) +H_key(3) + 1

          U2 = U(N2*2-1)*c +U(N2*2)*s
          U1 = U(N1*2-1)*c +U(N1*2)*s
          Def(2) = U2 - U1

      return
      end
          
      