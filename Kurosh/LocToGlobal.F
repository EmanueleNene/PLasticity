      subroutine LocToGlobal(ElNum,dXYZ_L,
     &                       ac,bc,dir,ta,XYZ_G,dXYZ_G)

#include "impcom.inc"
#include "beam_hyg.inc"
      integer             ElNum,nPln
      double precision    dXYZ_L(3),XYZ_G(3), dXYZ_G(3,3)
      double precision    ac(3),bc(3),dir(3),ta(3),tb(3)
      
      double precision    nxt(3)
      

      nPln = H_key(3)*(H_key(4)+H_key(9)+H_key(10))      
      if ( ElNum .le. nPln*H_key(12) ) then
          call v_move (ac(1),XYZ_G(1),3)
      else
          call v_move (bc(1),XYZ_G(1),3)         
      endif
      
      nxt(1) = dir(2) * ta(3) - dir(3) * ta(2)
      nxt(2) = dir(3) * ta(1) - dir(1) * ta(3)
      nxt(3) = dir(1) * ta(2) - dir(2) * ta(1)
      
      dXYZ_G(:,1) = dXYZ_L(1) * ta
      dXYZ_G(:,2) = dXYZ_L(2) * nxt
      dXYZ_G(:,3) = dXYZ_L(3) * dir
      
      return
      end