      subroutine FindEffDef(U,BeamProps,ac,bc,dir,ta,Ang,
     &                      Connectivity,BndEps,W)
c --- primary function: Finds strain at the middle of the bond 
c      TEMPORARY VERSION
      
#include "impcom.inc"
#include "beam_hyg.inc"

      double precision    U(H_Dof),BeamProps(11,2),Ang
      double precision    ac(3),bc(3),dir(3),ta(3),tb(3)
      integer             Connectivity(H_nElNode,H_nElem)
      
      double precision    BndEps(8,2),W(2)

      if (H_key(1) .eq. 2) then
          if (H_key(15) .eq. 1) then
              call FindDef2D_1(U,BeamProps,Ang,Connectivity,BndEps)
          elseif (H_key(15) .eq. 2) then
              call FindDef2D_2(U,BeamProps,Ang,Connectivity,BndEps,W)
          elseif (H_key(15) .eq. 3) then
              call FindDef2D_3(U,BeamProps,Ang,Connectivity,BndEps)
          endif

      elseif (H_key(1).eq. 3) then
          if (H_key(2) .eq. 1) then
              call FindDef3D_2(U,BeamProps,ac,bc,dir,ta,Ang,
     &                         Connectivity,BndEps)
c              call FindDef3D_2_1El(U,BeamProps,Ang,Connectivity,BndEps)
          elseif (H_key(2) .eq. 2) then
              call FindDef3D_2(U,BeamProps,ac,bc,dir,ta,Ang,
     &                         Connectivity,BndEps)
c              call FindDef3D_2_1El(U,BeamProps,Ang,Connectivity,BndEps)
          endif
      endif
            
      return
      end
          
      