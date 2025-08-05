      subroutine ElemLocX(ElNum,BeamProps,Ang,dXYZ)

#include "impcom.inc"
#include "beam_hyg.inc"

      integer             ElNum
      double precision    BeamProps(11,2), Ang
      double precision    dXYZ(3)
      integer             nPln 
      
      dXYZ(1) = BeamProps(3,2)/sin(Ang)/H_key(3)
      dXYZ(2) = BeamProps(3,1)/H_key(4)
      nPln = H_key(3)*(H_key(4)+H_key(9)+H_key(10))
      
      if (ElNum .le. nPln*H_key(12) ) then
          dXYZ(3) = BeamProps(2,1)/H_key(12)
      else
          dXYZ(3) = BeamProps(2,2)/H_key(11)
      endif
      
      dXYZ = 0.5d0 * dXYZ
      
      return
      end