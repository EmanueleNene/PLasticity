      subroutine FindDef3D1_1(U,Ang,Def)
c --- primary function: Finds deformation at the middle of the bond 

#include "impcom.inc"
#include "beam_hyg.inc"

      double precision    U(H_Dof),Ang
      double precision    Def(4)
      
      integer             N1, N2, N3, N4
      double precision    U1, U2
      double precision    c , s

      c = cos(Ang)
      s = sin(Ang)
      
      N1 = (H_key(3)+1) * H_key(9) + 
     &     (H_key(4)/2) *(H_key(7)+H_key(3)+H_key(8)+1) +
     &      H_key(7)+1
      N2 = N1 + H_key(3)
      if (2*(H_key(4)/2) .eq. H_key(4)) then
          Def(1) = U(N2*2-1) - U(N1*2-1)
      else
          N3 = N1 + H_key(7)+H_key(3)+H_key(8)+1
          N4 = N3 + H_key(3)
          Def(1) = 0.5d0* ( U(N2*2-1)+U(N4*2-1)
     &                     -U(N1*2-1)-U(N3*2-1) )
      endif
      N1 = (H_key(3)+1) * H_key(9) + 
     &      H_key(7) +H_key(3)/2 +1
      N2 = (H_key(3)+1) * H_key(9) + 
     &      H_key(4) *(H_key(7)+H_key(3)+H_key(8)+1) +
     &      H_key(7) +H_key(3)/2 +1
      if (2*(H_key(3)/2) .eq. H_key(3)) then
          U2 = U(N2*2-1)*c +U(N2*2)*s
          U1 = U(N1*2-1)*c +U(N1*2)*s
          Def(2) = U2 - U1
      else
          U2 = (U(N2*2-1)+U(N2*2+1))*c+(U(N2*2)+U(N2*2+2))*s
          U1 = (U(N1*2-1)+U(N1*2+1))*c+(U(N1*2)+U(N1*2+2))*s
          Def(2) = 0.5d0* ( U2 - U1 )
      endif
      
      
      
      Def(1) = 0.5d0*(U(16)-U(13) + U(22)-U(19));
      Def(3) = 0.5d0*(U(4)-U(1) + U(10)-U(7));
      
      Def(2) = 0.5d0*(U(16)-U(13) + U(22)-U(19)) * c +
     &         0.5d0*(U(20)-U(14) + U(23)-U(17)) * s
      Def(4) = 0.5d0*(U(31)-U(25) + U(34)-U(28)) * c +
     &         0.5d0*(U(29)-U(29) + U(35)-U(32)) * s
      

      return
      end
          
      