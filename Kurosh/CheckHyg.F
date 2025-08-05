      integer function CheckHyg(BeamProps,Ang)

#include "impcom.inc"
#include "solu_info.inc"
      
      double precision    BeamProps(9,2),Ang
      double precision    angtols,ONE
      double precision    BetaL1, BetaT1, DH1, BetaL2, BetaT2, DH2
      parameter          (angtols = 0.02d0,
     &                    ONE = 1.0d0)

      CheckHyg = 0
      
      BetaL1  = BeamProps(7,1)
      BetaT1  = BeamProps(8,1)
      DH1     = BeamProps(9,1)
      BetaL2  = BeamProps(7,2)
      BetaT2  = BeamProps(8,2)
      DH2     = BeamProps(9,2)
      
      if ( (BetaL1*DH1+ONE).le.TINY .or.
     &     (BetaT1*DH1+ONE).le.TINY .or.
     &     (BetaL2*DH2+ONE).le.TINY .or.
     &     (BetaT2*DH2+ONE).le.TINY ) then
          call erhandler ('BndHygStrn',5000,3,
     &        'Too much shrinkage, check hygroexpansion coefficient and 
     &         moisture change (%g). Hygro effects ignored!', 0.0d0,' ')
          CheckHyg = 2
          call StopAll(0)
          goto 999
      endif
      
      if (A_KeyFirstPass.eq.0) goto 999
      
      if (abs(sin(Ang)).lt.angtols) then
          call erhandler ('BndHygStrn',5000,2,
     &        'Contact angle is smaller than 2 degrees',0.0d0,' ')
          CheckHyg = 1
          goto 999
      endif


999   return
      end