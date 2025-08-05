      subroutine SetConnect3D1(Connectivity)
c --- primary function: setup the element connectivity for 2D linear elements
      
#include "impcom.inc"
#include "beam_hyg.inc"

      integer     Connectivity(H_nElNode,H_nElem)
      
      integer     i,j,k,count,c,cc,ccc,cccc,nLin,nPln

c     nLin: Number of nodes per line  as we proceed with meshing
c     nPln: Number of nodes per plane as we proceed with meshing
      
      c = 0
      count = 0
      nLin = H_key(3)+1                                   ! Number of nodes per line that we sweep (lines in x direction, sweep in y direction)
      nPln = (H_key(4)+H_key(9)+H_key(10) +1)*nLin        ! Number of nodes per plane that we sweep (in z direction)
      if (H_key(5) .eq. 0) then
          do i=1,H_key(12)
              do j=1,H_key(4)+H_key(9)+H_key(10)
                  do k=1,H_key(3)
                      count = count + 1
                      Connectivity(1,count)= (i-1)*nPln +(j-1)*nLin +k
                      Connectivity(2,count)= (i-1)*nPln +(j-1)*nLin +k+1
                      Connectivity(3,count)= (i-1)*nPln +(j  )*nLin +k+1
                      Connectivity(4,count)= (i-1)*nPln +(j  )*nLin +k
                      Connectivity(5,count)= (i  )*nPln +(j-1)*nLin +k
                      Connectivity(6,count)= (i  )*nPln +(j-1)*nLin +k+1
                      Connectivity(7,count)= (i  )*nPln +(j  )*nLin +k+1
                      Connectivity(8,count)= (i  )*nPln +(j  )*nLin +k
                    enddo
              enddo
          enddo
          c= nPln*(H_key(12)+1)
      elseif (H_key(14) .eq. 0) then
      else
      endif
      
      nLin = H_key(4)+1
      nPln = (H_key(3)+H_key(7)+H_key(8) +1)*nLin
      cc   = c - (H_key(3)+1)*(H_key(4)+1)
      ccc  = (H_key(7)-1) * nLin
      cccc = cc - (H_key(3)+1)*H_key(10) 
      if (H_key(7) .gt. 0) then
          do j=1,H_key(7)-1
            do k=1,H_key(4)
                count = count + 1
                Connectivity(1,count)= c+       +(j-1)*nLin +k
                Connectivity(2,count)= c+       +(j  )*nLin +k
                Connectivity(3,count)= c+       +(j  )*nLin +k+1
                Connectivity(4,count)= c+       +(j-1)*nLin +k+1
                Connectivity(5,count)= cc+ nPln +(j-1)*nLin +k
                Connectivity(6,count)= cc+ nPln +(j  )*nLin +k
                Connectivity(7,count)= cc+ nPln +(j  )*nLin +k+1
                Connectivity(8,count)= cc+ nPln +(j-1)*nLin +k+1
            enddo
          enddo
          do k=1,H_key(4)
              count = count + 1
              Connectivity(1,count)= c+       +ccc        +k
              Connectivity(2,count)= cccc + (k-1)*(H_key(3)+1) +1
              Connectivity(3,count)= cccc + (k  )*(H_key(3)+1) +1
              Connectivity(4,count)= c+       +ccc        +k+1
              Connectivity(5,count)= cc+ nPln +ccc +      +k
              Connectivity(6,count)= cc+ nPln +ccc + nLin +k
              Connectivity(7,count)= cc+ nPln +ccc + nLin +k+1
              Connectivity(8,count)= cc+ nPln +ccc +      +k+1
          enddo
      endif
      cccc = cccc - H_key(7)
      do j=H_key(7)+1,H_key(7)+H_key(3)
        do k=1,H_key(4)
            count = count + 1
            Connectivity(1,count)= cccc + (k-1)*(H_key(3)+1) + j
            Connectivity(2,count)= cccc + (k-1)*(H_key(3)+1) + j+1
            Connectivity(3,count)= cccc + (k  )*(H_key(3)+1) + j+1
            Connectivity(4,count)= cccc + (k  )*(H_key(3)+1) + j
            Connectivity(5,count)= cc+ nPln +(j-1)*nLin +k
            Connectivity(6,count)= cc+ nPln +(j  )*nLin +k
            Connectivity(7,count)= cc+ nPln +(j  )*nLin +k+1
            Connectivity(8,count)= cc+ nPln +(j-1)*nLin +k+1
        enddo
      enddo
      if (H_key(8) .gt. 0) then
        cccc = cccc+H_key(7)+H_key(3)
        ccc  = ccc + (H_key(3)+1)* nLin
        do k=1,H_key(4)
            count = count + 1
            Connectivity(1,count)= cccc + (k-1)*(H_key(3)+1) +1
            Connectivity(2,count)= cc+      +(j  )*nLin +k
            Connectivity(3,count)= cc+      +(j  )*nLin +k+1
            Connectivity(4,count)= cccc + (k  )*(H_key(3)+1) +1
            Connectivity(5,count)= cc+ nPln +ccc      +k
            Connectivity(6,count)= cc+ nPln +ccc+nLin +k
            Connectivity(7,count)= cc+ nPln +ccc+nLin +k+1
            Connectivity(8,count)= cc+ nPln +ccc      +k+1
        enddo
        do j=H_key(7)+H_key(3)+2,H_key(7)+H_key(3)+H_key(8)
            do k=1,H_key(4)
                count = count + 1
                Connectivity(1,count)= cc+      +(j-1)*nLin +k
                Connectivity(2,count)= cc+      +(j  )*nLin +k
                Connectivity(3,count)= cc+      +(j  )*nLin +k+1
                Connectivity(4,count)= cc+      +(j-1)*nLin +k+1
                Connectivity(5,count)= cc+ nPln +(j-1)*nLin +k
                Connectivity(6,count)= cc+ nPln +(j  )*nLin +k
                Connectivity(7,count)= cc+ nPln +(j  )*nLin +k+1
                Connectivity(8,count)= cc+ nPln +(j-1)*nLin +k+1
            enddo
          enddo
      endif
      
      if (H_key(6) .eq. 0) then
          do i=2,H_key(11)
              do j=1,H_key(3)+H_key(7)+H_key(8)
                  do k=1,H_key(4)
                      count = count + 1
                      Connectivity(1,count)=cc+(i-1)*nPln+(j-1)*nLin+k
                      Connectivity(2,count)=cc+(i-1)*nPln+(j  )*nLin+k
                      Connectivity(3,count)=cc+(i-1)*nPln+(j  )*nLin+k+1
                      Connectivity(4,count)=cc+(i-1)*nPln+(j-1)*nLin+k+1
                      Connectivity(5,count)=cc+(i  )*nPln+(j-1)*nLin+k
                      Connectivity(6,count)=cc+(i  )*nPln+(j  )*nLin+k
                      Connectivity(7,count)=cc+(i  )*nPln+(j  )*nLin+k+1
                      Connectivity(8,count)=cc+(i  )*nPln+(j-1)*nLin+k+1
                    enddo
              enddo
          enddo
      elseif (H_key(14) .eq. 0) then
      else
      endif
      
999   return
      end