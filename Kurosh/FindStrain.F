      subroutine FindStrain(BeamProps,Xi,Def,Eps)
#include "impcom.inc"

      
      double precision    BeamProps(9,2), Xi(2), Def(2)
      double precision    Eps(2,2)
      
      double precision    L1,L2
      
      
      L1      = BeamProps(1,1)
      L2      = BeamProps(1,2)
      
      Eps(1,1) = Def(1)/L1*(1.0d0-3.0d0*Xi(1))
      Eps(2,1) = Def(1)/L1*(1.0d0+3.0d0*Xi(1))
      
      Eps(1,2) = Def(2)/L2*(1.0d0-3.0d0*Xi(2))
      Eps(2,2) = Def(2)/L2*(1.0d0+3.0d0*Xi(2))
      
      return
      end