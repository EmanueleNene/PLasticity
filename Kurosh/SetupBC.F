      subroutine SetupBC()
c --- primary function: Sets the boundary conditions
      
#include "impcom.inc"
#include "beam_hyg.inc"
#include "constants.inc"

      PTRFTN  MemAlloc
      integer             i,j,count,count0
      
      H_pBCData = MemAlloc(H_Dof*4,W_LenInt,'H_pBCData       ')
      
      if (H_key(1) .eq. 2) then
          H_nBC = 3
          if (H_key(2) .eq. 1) then
              H_BCData(1) = (H_key(7) +  H_key(9)*   (H_key(3)+1)) *2+ 1
              H_BCData(2) = (H_key(7) +  H_key(9)*   (H_key(3)+1)) *2+ 2
              H_BCData(3) = (H_key(7) + (H_key(9)+1)*(H_key(3)+1)) *2
          elseif (H_key(2) .eq. 2) then
              H_BCData(1) = ( H_key(9)*(3*H_key(3)+2) +2*H_key(7) ) *2+1
              H_BCData(2) = ( H_key(9)*(3*H_key(3)+2) +2*H_key(7) ) *2+2
              H_BCData(3) = ( H_key(9)*(3*H_key(3)+2) +
     &                                   2*(H_key(3)+H_key(7)) ) *2 +2
          endif
      elseif (H_key(1) .eq. 3) then
          H_nBC = 6
          if (H_key(2) .eq. 1) then
              i = (H_key(3)+1)*(H_key(4)+H_key(9)+H_key(10)+1)*H_key(12)
     &           +(H_key(3)+1) * H_key(9) + 1
              H_BCData(1) = i*3-2
              H_BCData(2) = i*3-1
              H_BCData(3) = i*3
              j = i + H_key(3)
              H_BCData(4) = j*3-1
              H_BCData(5) = j*3
              j = i + (H_key(3)+1)*H_key(4)
              H_BCData(6) = j*3
          elseif (H_key(2) .eq. 2) then
              i = ( (4*H_key(3)+3) * (H_key(4)+H_key(9)+H_key(10)) +
     &              (3*H_key(3)+2) ) *H_key(12)
     &           +(3*H_key(3)+2)*H_key(9) + 1
              H_BCData(1) = i*3-2
              H_BCData(2) = i*3-1
              H_BCData(3) = i*3
              j = i + 2*H_key(3)
              H_BCData(4) = j*3-1
              H_BCData(5) = j*3
              j = i + (3*H_key(3)+2)*H_key(4)
              H_BCData(6) = j*3
          endif
      endif
      
      H_pBCLoc   = 1
      H_pFreeLoc = H_pBCLoc + H_nBC
      
      count0 = 0
      do i = 1, H_nBC
          count = count0
          do j=count0+i,H_BCData(i)-1
              count = count+1
              H_BCData(count+H_nBC) = j
          enddo
          count0 = count
      enddo
      do j=count0+i,H_Dof
          count = count+1
          H_BCData(count+H_nBC) = j
      enddo
      
      return
      end
