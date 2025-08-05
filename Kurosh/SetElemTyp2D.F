      subroutine SetElemTyp2D(ElemType)
c --- primary function: setup the element stiffness types for 2D elements
      
#include "impcom.inc"
#include "beam_hyg.inc"
      
      integer     ElemType(H_nElem)
      
      integer     c,count,i,j
      
      c = 0
      count = 0
      do i=1,H_key(9)
          do j=1,H_key(3)
              count = count + 1
                            
              if ( (j .gt. (H_key(3)-H_key(5))/2) .and.
     &             (j .le. (H_key(3)+H_key(5))/2) )then
                ElemType(count) = 14
              else
                  ElemType(count) = 13
              endif

          enddo
      enddo
                  
      do i=1,H_key(4)
          do j=1,H_key(7)+H_key(3)+H_key(8)
              count = count + 1
                          
              if  (j .le. H_key(7)) then
                  if ( (i .gt. (H_key(4)-H_key(6))/2) .and.
     &                 (i .le. (H_key(4)+H_key(6))/2) ) then
                      ElemType(count) = 10
                  else
                      ElemType(count) = 9
                  endif
              elseif (j .gt. H_key(3)+H_key(7)) then
                  if ( (i .gt. (H_key(4)-H_key(6))/2) .and.
     &                 (i .le. (H_key(4)+H_key(6))/2) ) then
                      ElemType(count) = 12
                  else
                      ElemType(count) = 11
                  endif
              elseif ((j.gt.(H_key(7)+(H_key(3)-H_key(5))/2)).and.
     &                (j.le.(H_key(7)+(H_key(3)+H_key(5))/2)))then
                  if ( (i .gt. (H_key(4)-H_key(6))/2) .and.
     &                 (i .le. (H_key(4)+H_key(6))/2) ) then
                      ElemType(count) = 7
                      count = count + 1
                      ElemType(count) = 8
                  else
                      ElemType(count) = 5
                      count = count + 1
                      ElemType(count) = 6
                  endif
              else
                  if ( (i .gt. (H_key(4)-H_key(6))/2) .and.
     &                 (i .le. (H_key(4)+H_key(6))/2) ) then
                      ElemType(count) = 3
                      count = count + 1
                      ElemType(count) = 4
                  else
                      ElemType(count) = 1
                      count = count + 1
                      ElemType(count) = 2
                  endif
              endif
                          
          enddo
      enddo
                  
      do i=1,H_key(10)
          do j=1,H_key(3)
              count = count + 1
              
              if ( (j .gt. (H_key(3)-H_key(5))/2) .and.
     &             (j .le. (H_key(3)+H_key(5))/2) )then
                ElemType(count) = 16
              else
                  ElemType(count) = 15
              endif
              
          enddo
      enddo

      return
      end