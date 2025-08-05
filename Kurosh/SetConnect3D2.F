      subroutine SetConnect3D2(Connectivity)
c --- primary function: setup the element connectivity for 2D linear elements
      
#include "impcom.inc"
#include "beam_hyg.inc"

      integer     Connectivity(H_nElNode,H_nElem)
      
      integer     i,j,k,count
      integer     c,c2,c3,cc,ccc,nComNod
      integer     nLin,n2Lin,nPln,n2Pln,nLin0,n2Lin0

c     nLin: Number of nodes per line  as we proceed with meshing
c     nPln: Number of nodes per plane as we proceed with meshing
      
      c = 0
      count = 0
      nLin  = 2*H_key(3)+1                                            ! Number of nodes per edge line that we sweep (lines in x direction, sweep in y direction)
      n2Lin = nLin + H_key(3)+1                                       ! Number of nodes per one edge and one mid line that we sweep (lines in x direction, sweep in y direction)
      nPln  = (H_key(4)+H_key(9)+H_key(10))*n2Lin + nLin              ! Number of nodes per edge plane that we sweep (in z direction)
      n2Pln = nPln + (H_key(4)+H_key(9)+H_key(10)+1)*(H_key(3)+1)     ! Number of nodes per edge plane that we sweep (in z direction)
      if (H_key(5) .eq. 0) then
          do i=1,H_key(12)
              do j=1,H_key(4)+H_key(9)+H_key(10)
                  do k=1,H_key(3)
                    count = count + 1
       Connectivity(1 ,count)= (i-1)*n2Pln +(j-1)*n2Lin +(k-1)*2+1
       Connectivity(2 ,count)= (i-1)*n2Pln +(j-1)*n2Lin +(k  )*2+1
       Connectivity(3 ,count)= (i-1)*n2Pln +(j  )*n2Lin +(k  )*2+1
       Connectivity(4 ,count)= (i-1)*n2Pln +(j  )*n2Lin +(k-1)*2+1
       Connectivity(5 ,count)= (i  )*n2Pln +(j-1)*n2Lin +(k-1)*2+1
       Connectivity(6 ,count)= (i  )*n2Pln +(j-1)*n2Lin +(k  )*2+1
       Connectivity(7 ,count)= (i  )*n2Pln +(j  )*n2Lin +(k  )*2+1
       Connectivity(8 ,count)= (i  )*n2Pln +(j  )*n2Lin +(k-1)*2+1
                    
       Connectivity(9 ,count)= (i-1)*n2Pln +(j-1)*n2Lin +k+nLin+1
       Connectivity(10,count)= (i-1)*n2Pln +(j  )*n2Lin +(k  )*2
       Connectivity(11,count)= (i-1)*n2Pln +(j-1)*n2Lin +k+nLin
       Connectivity(12,count)= (i-1)*n2Pln +(j-1)*n2Lin +(k  )*2
       Connectivity(13,count)= (i  )*n2Pln +(j-1)*n2Lin +k+nLin+1
       Connectivity(14,count)= (i  )*n2Pln +(j  )*n2Lin +(k  )*2
       Connectivity(15,count)= (i  )*n2Pln +(j-1)*n2Lin +k+nLin
       Connectivity(16,count)= (i  )*n2Pln +(j-1)*n2Lin +(k  )*2
            
       Connectivity(17,count)= (i-1)*n2Pln+nPln +(j-1)*(n2Lin-nLin) +k+1
       Connectivity(18,count)= (i-1)*n2Pln+nPln +(j  )*(n2Lin-nLin) +k+1
       Connectivity(19,count)= (i-1)*n2Pln+nPln +(j  )*(n2Lin-nLin) +k
       Connectivity(20,count)= (i-1)*n2Pln+nPln +(j-1)*(n2Lin-nLin) +k
                    
                    enddo
              enddo
          enddo
          c  = n2Pln*H_key(12)+nPln
          c2 = n2Pln*H_key(12) + n2Lin*H_key(9)
      elseif (H_key(14) .eq. 0) then
      else
      endif
      
      nLin0  = nLin
      n2Lin0 = n2Lin
      nLin  = 2*H_key(4)+1
      n2Lin = nLin + H_key(4)+1
      nPln  = (H_key(3)+H_key(7)+H_key(8))*n2Lin + nLin
      n2Pln = nPln + (H_key(3)+H_key(7)+H_key(8)+1)*(H_key(4)+1)
      
      nComNod = (2*H_key(3)+1)*(2*H_key(4)+1)-H_key(3)*H_key(4)
      cc      = c + n2Pln - nComNod
      ccc     = cc - (n2Pln-nPln)
      if (H_key(7) .gt. 0) then
          do j=1,H_key(7)-1
            do k=1,H_key(4)
                count = count + 1
                Connectivity(1 ,count)= c +  (j-1)*n2Lin +(k-1)*2+1
                Connectivity(2 ,count)= c +  (j  )*n2Lin +(k-1)*2+1
                Connectivity(3 ,count)= c +  (j  )*n2Lin +(k  )*2+1
                Connectivity(4 ,count)= c +  (j-1)*n2Lin +(k  )*2+1
                Connectivity(5 ,count)= cc+  (j-1)*n2Lin +(k-1)*2+1
                Connectivity(6 ,count)= cc+  (j  )*n2Lin +(k-1)*2+1
                Connectivity(7 ,count)= cc+  (j  )*n2Lin +(k  )*2+1
                Connectivity(8 ,count)= cc+  (j-1)*n2Lin +(k  )*2+1
                
                Connectivity(9 ,count)= c +  (j  )*n2Lin +(k  )*2
                Connectivity(10,count)= c +  (j-1)*n2Lin +k+nLin+1
                Connectivity(11,count)= c +  (j-1)*n2Lin +(k  )*2
                Connectivity(12,count)= c +  (j-1)*n2Lin +k+nLin
                Connectivity(13,count)= cc+  (j  )*n2Lin +(k  )*2
                Connectivity(14,count)= cc+  (j-1)*n2Lin +k+nLin+1
                Connectivity(15,count)= cc+  (j-1)*n2Lin +(k  )*2
                Connectivity(16,count)= cc+  (j-1)*n2Lin +k+nLin
                
                Connectivity(17,count)= ccc+ (j  )*(n2Lin-nLin) +k
                Connectivity(18,count)= ccc+ (j  )*(n2Lin-nLin) +k+1
                Connectivity(19,count)= ccc+ (j-1)*(n2Lin-nLin) +k+1
                Connectivity(20,count)= ccc+ (j-1)*(n2Lin-nLin) +k
            enddo
          enddo
          
          j    = H_key(7)
          do k=1,H_key(4)
              count = count + 1
              Connectivity(1 ,count)= c +  (j-1)*n2Lin  +(k-1)*2+1
              Connectivity(2 ,count)= c2+  (k-1)*n2Lin0 +1
              Connectivity(3 ,count)= c2+  (k  )*n2Lin0 +1
              Connectivity(4 ,count)= c +  (j-1)*n2Lin +(k  )*2+1
              Connectivity(5 ,count)= cc+  (j-1)*n2Lin +(k-1)*2+1
              Connectivity(6 ,count)= cc+  (j  )*n2Lin +(k-1)*2+1
              Connectivity(7 ,count)= cc+  (j  )*n2Lin +(k  )*2+1
              Connectivity(8 ,count)= cc+  (j-1)*n2Lin +(k  )*2+1
               
              Connectivity(9 ,count)= c2+  (k-1)*n2Lin0 +nLin0+1
              Connectivity(10,count)= c +  (j-1)*n2Lin +k+nLin+1
              Connectivity(11,count)= c +  (j-1)*n2Lin +(k  )*2
              Connectivity(12,count)= c +  (j-1)*n2Lin +k+nLin
              Connectivity(13,count)= cc+  (j  )*n2Lin +(k  )*2
              Connectivity(14,count)= cc+  (j-1)*n2Lin +k+nLin+1
              Connectivity(15,count)= cc+  (j-1)*n2Lin +(k  )*2
              Connectivity(16,count)= cc+  (j-1)*n2Lin +k+nLin
                
              Connectivity(17,count)= ccc+ (j  )*(n2Lin-nLin) +k
              Connectivity(18,count)= ccc+ (j  )*(n2Lin-nLin) +k+1
              Connectivity(19,count)= ccc+ (j-1)*(n2Lin-nLin) +k+1
              Connectivity(20,count)= ccc+ (j-1)*(n2Lin-nLin) +k
          enddo
      endif
      
      do j=H_key(7)+1,H_key(7)+H_key(3)
        i = j-H_key(7)
        do k=1,H_key(4)
            count = count + 1
          Connectivity(1 ,count)= c2+  (k-1)*n2Lin0 + (i-1)*2+1
          Connectivity(2 ,count)= c2+  (k-1)*n2Lin0 + (i  )*2+1
          Connectivity(3 ,count)= c2+  (k  )*n2Lin0 + (i  )*2+1
          Connectivity(4 ,count)= c2+  (k  )*n2Lin0 + (i-1)*2+1
          Connectivity(5 ,count)= cc+  (j-1)*n2Lin +(k-1)*2+1
          Connectivity(6 ,count)= cc+  (j  )*n2Lin +(k-1)*2+1
          Connectivity(7 ,count)= cc+  (j  )*n2Lin +(k  )*2+1
          Connectivity(8 ,count)= cc+  (j-1)*n2Lin +(k  )*2+1
               
          Connectivity(9 ,count)= c2+  (k-1)*n2Lin0 +i+nLin0+1
          Connectivity(10,count)= c2+  (k  )*n2Lin0 +i*2
          Connectivity(11,count)= c2+  (k-1)*n2Lin0 +i+nLin0
          Connectivity(12,count)= c2+  (k-1)*n2Lin0 +i*2
          Connectivity(13,count)= cc+  (j  )*n2Lin +(k  )*2
          Connectivity(14,count)= cc+  (j-1)*n2Lin +k+nLin+1
          Connectivity(15,count)= cc+  (j-1)*n2Lin +(k  )*2
          Connectivity(16,count)= cc+  (j-1)*n2Lin +k+nLin
                
          Connectivity(17,count)= ccc+ (j  )*(n2Lin-nLin) +k
          Connectivity(18,count)= ccc+ (j  )*(n2Lin-nLin) +k+1
          Connectivity(19,count)= ccc+ (j-1)*(n2Lin-nLin) +k+1
          Connectivity(20,count)= ccc+ (j-1)*(n2Lin-nLin) +k
        enddo
      enddo
      if (H_key(8) .gt. 0) then
        c3 = c -nComNod
        i    = H_key(3)
        do k=1,H_key(4)
            count = count + 1
            Connectivity(1 ,count)= c2+  (k-1)*n2Lin0 + (i  )*2+1
            Connectivity(2 ,count)= c3 +  (j  )*n2Lin +(k-1)*2+1
            Connectivity(3 ,count)= c3 +  (j  )*n2Lin +(k  )*2+1
            Connectivity(4 ,count)= c2+  (k  )*n2Lin0 + (i  )*2+1
            Connectivity(5 ,count)= cc+  (j-1)*n2Lin +(k-1)*2+1
            Connectivity(6 ,count)= cc+  (j  )*n2Lin +(k-1)*2+1
            Connectivity(7 ,count)= cc+  (j  )*n2Lin +(k  )*2+1
            Connectivity(8 ,count)= cc+  (j-1)*n2Lin +(k  )*2+1
            
            Connectivity(9 ,count)= c3+  (j  )*n2Lin +(k  )*2
            Connectivity(10,count)= c3+  (j-1)*n2Lin +k+nLin+1
            Connectivity(11,count)= c2+  (k-1)*n2Lin0 +n2Lin0
            Connectivity(12,count)= c3+  (j-1)*n2Lin +k+nLin
            
            Connectivity(13,count)= cc+  (j  )*n2Lin +(k  )*2
            Connectivity(14,count)= cc+  (j-1)*n2Lin +k+nLin+1
            Connectivity(15,count)= cc+  (j-1)*n2Lin +(k  )*2
            Connectivity(16,count)= cc+  (j-1)*n2Lin +k+nLin
                
            Connectivity(17,count)= ccc+ (j  )*(n2Lin-nLin) +k
            Connectivity(18,count)= ccc+ (j  )*(n2Lin-nLin) +k+1
            Connectivity(19,count)= ccc+ (j-1)*(n2Lin-nLin) +k+1
            Connectivity(20,count)= ccc+ (j-1)*(n2Lin-nLin) +k
        enddo
        do j=H_key(7)+H_key(3)+2,H_key(7)+H_key(3)+H_key(8)
            do k=1,H_key(4)
                count = count + 1
                Connectivity(1 ,count)= c3+  (j-1)*n2Lin +(k-1)*2+1
                Connectivity(2 ,count)= c3+  (j  )*n2Lin +(k-1)*2+1
                Connectivity(3 ,count)= c3+  (j  )*n2Lin +(k  )*2+1
                Connectivity(4 ,count)= c3+  (j-1)*n2Lin +(k  )*2+1
                Connectivity(5 ,count)= cc+  (j-1)*n2Lin +(k-1)*2+1
                Connectivity(6 ,count)= cc+  (j  )*n2Lin +(k-1)*2+1
                Connectivity(7 ,count)= cc+  (j  )*n2Lin +(k  )*2+1
                Connectivity(8 ,count)= cc+  (j-1)*n2Lin +(k  )*2+1
                
                Connectivity(9 ,count)= c3+  (j  )*n2Lin +(k  )*2
                Connectivity(10,count)= c3+  (j-1)*n2Lin +k+nLin+1
                Connectivity(11,count)= c3+  (j-1)*n2Lin +(k  )*2
                Connectivity(12,count)= c3+  (j-1)*n2Lin +k+nLin
                Connectivity(13,count)= cc+  (j  )*n2Lin +(k  )*2
                Connectivity(14,count)= cc+  (j-1)*n2Lin +k+nLin+1
                Connectivity(15,count)= cc+  (j-1)*n2Lin +(k  )*2
                Connectivity(16,count)= cc+  (j-1)*n2Lin +k+nLin
                
                Connectivity(17,count)= ccc+ (j  )*(n2Lin-nLin) +k
                Connectivity(18,count)= ccc+ (j  )*(n2Lin-nLin) +k+1
                Connectivity(19,count)= ccc+ (j-1)*(n2Lin-nLin) +k+1
                Connectivity(20,count)= ccc+ (j-1)*(n2Lin-nLin) +k
            enddo
          enddo
      endif
      
      if (H_key(6) .eq. 0) then
          do i=1,H_key(11)-1
              do j=1,H_key(3)+H_key(7)+H_key(8)
                  do k=1,H_key(4)
                      count = count + 1
       Connectivity(1 ,count)= cc +(i-1)*n2Pln +(j-1)*n2Lin +(k-1)*2+1
       Connectivity(2 ,count)= cc +(i-1)*n2Pln +(j  )*n2Lin +(k-1)*2+1
       Connectivity(3 ,count)= cc +(i-1)*n2Pln +(j  )*n2Lin +(k  )*2+1
       Connectivity(4 ,count)= cc +(i-1)*n2Pln +(j-1)*n2Lin +(k  )*2+1
       Connectivity(5 ,count)= cc +(i  )*n2Pln +(j-1)*n2Lin +(k-1)*2+1
       Connectivity(6 ,count)= cc +(i  )*n2Pln +(j  )*n2Lin +(k-1)*2+1
       Connectivity(7 ,count)= cc +(i  )*n2Pln +(j  )*n2Lin +(k  )*2+1
       Connectivity(8 ,count)= cc +(i  )*n2Pln +(j-1)*n2Lin +(k  )*2+1
                    
       Connectivity(9 ,count)= cc +(i-1)*n2Pln +(j  )*n2Lin +(k  )*2
       Connectivity(10,count)= cc +(i-1)*n2Pln +(j-1)*n2Lin +k+nLin+1
       Connectivity(11,count)= cc +(i-1)*n2Pln +(j-1)*n2Lin +(k  )*2
       Connectivity(12,count)= cc +(i-1)*n2Pln +(j-1)*n2Lin +k+nLin
       Connectivity(13,count)= cc +(i  )*n2Pln +(j  )*n2Lin +(k  )*2
       Connectivity(14,count)= cc +(i  )*n2Pln +(j-1)*n2Lin +k+nLin+1
       Connectivity(15,count)= cc +(i  )*n2Pln +(j-1)*n2Lin +(k  )*2
       Connectivity(16,count)= cc +(i  )*n2Pln +(j-1)*n2Lin +k+nLin
       
            
      Connectivity(17,count)=cc+(i-1)*n2Pln+nPln +(j  )*(n2Lin-nLin)+k
      Connectivity(18,count)=cc+(i-1)*n2Pln+nPln +(j  )*(n2Lin-nLin)+k+1
      Connectivity(19,count)=cc+(i-1)*n2Pln+nPln +(j-1)*(n2Lin-nLin)+k+1
      Connectivity(20,count)=cc+(i-1)*n2Pln+nPln +(j-1)*(n2Lin-nLin)+k
                    enddo
              enddo
          enddo
      elseif (H_key(14) .eq. 0) then
      else
      endif
      
999   return
      end