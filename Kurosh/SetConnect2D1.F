      subroutine SetConnect2D1(Connectivity)
c --- primary function: setup the element connectivity for 2D linear elements
      
#include "impcom.inc"
#include "beam_hyg.inc"

      integer     Connectivity(H_nElNode,H_nElem)
      
      integer     i,j,count,c,cc,nVerNod
      
      c = 0
      count = 0
      nVerNod = H_key(3)+1
      do i=1,H_key(9)
          do j=1,H_key(3)
              count = count + 1
              if (i .eq. H_key(9)) then
                  c = H_key(7)
              endif
                          
              Connectivity(1,count) = (i-1)*(nVerNod) + j
              Connectivity(2,count) = (i-1)*(nVerNod) + j+1
              Connectivity(3,count) = c + i*(nVerNod) + j+1
              Connectivity(4,count) = c + i*(nVerNod) + j
          enddo
      enddo
                  
      c = (H_key(3)+1) * H_key(9)
      nVerNod = H_key(7)+H_key(3)+H_key(8)+1
      do i=1,H_key(4)
          do j=1,H_key(7)+H_key(3)+H_key(8)
              count = count + 1
                          
              Connectivity(1,count) = c+(i-1)*(nVerNod)+j
              Connectivity(2,count) = c+(i-1)*(nVerNod)+j+1
              Connectivity(3,count) = c+    i*(nVerNod)+j+1
              Connectivity(4,count) = c+    i*(nVerNod)+j

              if ( (j .gt. H_key(7)) .and.
     &             (j .le. H_key(3)+H_key(7)) ) then
                  Connectivity(:,count+1) = Connectivity(:,count)
                  count = count + 1
              endif

          enddo
      enddo
                  
      c = (H_key(3)+1)*H_key(9) -H_key(3)-1 +
     &    (H_key(3)+H_key(7)+H_key(8)+1)*(H_key(4)+1)
      nVerNod = H_key(3)+1
      do i=1,H_key(10)
          do j=1,H_key(3)
              count = count + 1
              if (i .eq. 1) then
                  cc = c-H_key(8)
              else
                  cc = c
              endif
                          
              Connectivity(1,count) = cc+(i-1)*(nVerNod)+j
              Connectivity(2,count) = cc+(i-1)*(nVerNod)+j+1
              Connectivity(3,count) = c +    i*(nVerNod)+j+1
              Connectivity(4,count) = c +    i*(nVerNod)+j
          enddo
      enddo
      
999   return
      end