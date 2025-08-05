      subroutine SetConnect2D2(Connectivity)
c --- primary function: setup the element connectivity for 2D linear elements
      
#include "impcom.inc"
#include "beam_hyg.inc"

      integer     Connectivity(H_nElNode,H_nElem)
      
      integer     i,j,count,c,cc,nNodLay
      
      c = 0
      count = 0
      nNodLay = 3*H_key(3)+2
      do i=1,H_key(9)
          do j=1,H_key(3)
              count = count + 1
              if (i .eq. H_key(9)) then
                  c = 2*H_key(7)
              endif
                          
              Connectivity(1,count) = (i-1)*(nNodLay) + 2*j -1
              Connectivity(2,count) = (i-1)*(nNodLay) + 2*j +1
              Connectivity(3,count) =  i   *(nNodLay) + 2*j +1 +c
              Connectivity(4,count) =  i   *(nNodLay) + 2*j -1 +c
              Connectivity(5,count) = (i-1)*(nNodLay) + 2*j
              Connectivity(6,count) =  i   *(nNodLay) +   j -H_key(3)
              Connectivity(7,count) =  i   *(nNodLay) + 2*j    +c
              Connectivity(8,count) =  i   *(nNodLay) +   j -H_key(3)-1
          enddo
      enddo
                  
      c = (3*H_key(3)+2) * H_key(9)
      cc = H_key(7)+H_key(3)+H_key(8)
      nNodLay = 3*cc+2
      do i=1,H_key(4)
          do j=1,H_key(7)+H_key(3)+H_key(8)
              count = count + 1
              
              Connectivity(1,count) = c+ (i-1)*(nNodLay) + 2*j -1
              Connectivity(2,count) = c+ (i-1)*(nNodLay) + 2*j +1
              Connectivity(3,count) = c+  i   *(nNodLay) + 2*j +1
              Connectivity(4,count) = c+  i   *(nNodLay) + 2*j -1
              Connectivity(5,count) = c+ (i-1)*(nNodLay) + 2*j
              Connectivity(6,count) = c+  i   *(nNodLay) +   j -cc
              Connectivity(7,count) = c+  i   *(nNodLay) + 2*j
              Connectivity(8,count) = c+  i   *(nNodLay) +   j -cc-1
              
              if ( (j .gt. H_key(7)) .and.
     &             (j .le. H_key(3)+H_key(7)) ) then
                  Connectivity(:,count+1) = Connectivity(:,count)
                  count = count + 1
              endif

          enddo
      enddo
                  
      c = (3*H_key(3)+2)*(H_key(9)-1) +2*H_key(4) +1 +
     &    (H_key(3)+H_key(7)+H_key(8))*(3*H_key(4)+2) +H_key(3)+1
      nNodLay = 3*H_key(3)+2
      do i=1,H_key(10)
          do j=1,H_key(3)
              count = count + 1
              if (i .eq. 1) then
                  cc = c-2*H_key(8)
              else
                  cc = c
              endif
              
            Connectivity(1,count) = cc+ (i-1)*(nNodLay) + 2*j -1
            Connectivity(2,count) = cc+ (i-1)*(nNodLay) + 2*j +1
            Connectivity(3,count) = c +  i   *(nNodLay) + 2*j +1
            Connectivity(4,count) = c +  i   *(nNodLay) + 2*j -1
            Connectivity(5,count) = cc+ (i-1)*(nNodLay) + 2*j
            Connectivity(6,count) = c +  i   *(nNodLay) +   j-H_key(3)
            Connectivity(7,count) = c +  i   *(nNodLay) + 2*j
            Connectivity(8,count) = c +  i   *(nNodLay) +   j-H_key(3)-1
          enddo
      enddo
      
999   return
      end