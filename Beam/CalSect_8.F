      subroutine CalSect_8 (nnodi,numnodes,numele,elemconn,
     &                      len_stiff,address,column_heights,
     &                      nSecRstNodes)
c --- primary function: calculate section data
c
#include "impcom.inc"
#include "constants.inc"
c
      integer nnodi,MNODI
      parameter (MNODI = 9)
      integer numele, len_stiff, numnodes,
     &        column_heights(numnodes),
     &        address(numnodes),elemconn(MNODI,numele),
     &        j, i, k, n, l, jel, ls, nSecRstNodes

c --- assume one DOF per node in calculating the skyline
c
      call vi_zero(address(1),numnodes)

c --- the nSecRstNodes calculation here is for storing the number of
c     nodes where results will be stored
c
      nSecRstNodes = 0
      do k = 1,numele
         do l = 1,4
            address(elemconn(l,k)) = 1
         enddo
      enddo
      do k=1,numnodes
         if (address(k).eq.1) nSecRstNodes=nSecRstNodes+1
      enddo
      call vi_zero(address(1),numnodes)
      column_heights(1) = 1
      do i = 2,numnodes
         ls = 100000
         do k=1,numele
            jel = 0
            do l=1,nnodi
               if ((elemconn(l,k)- i).eq.0) then
                  jel = 1
               endif
            enddo
            if (jel.eq.1) then
                do n=1,nnodi
                   if ((elemconn(n,k)-ls) .lt. 0) then
                       ls = elemconn(n,k)
                   endif
                enddo
            endif
         enddo
         if (ls.eq.100000) ls = 0
         column_heights(i) = ls
      enddo
      do i=1,numnodes
         column_heights(i) = i-column_heights(i)
      enddo
      address(1) = 1
      len_stiff = column_heights(1)
      do j = 1,numnodes
         len_stiff = len_stiff+column_heights(j)
         address(j) = len_stiff+j
      enddo
      len_stiff = address(numnodes)
      
      return
      end
