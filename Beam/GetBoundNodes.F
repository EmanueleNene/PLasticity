      subroutine GetBoundNodes(ElemConn, numNodes, numElems, bNodes)
c --- primary function: given a cross section connectivity, return the boundary nodes
c
#include "impcom.inc"
c
      integer ElemConn(9,*),numElems, numNodes, bNodes(*), 
     &        i, j, k, l, m, n
c

#if defined(DEBUG)
      call SubBeg('GetBoundNodes')
#endif

c --- zero everything
      do i = 1,numNodes
         bNodes(i) = 0
      enddo
c
c --- mark all midside nodes
      do i = 1,numElems
         do j = 5,8
            k = ElemConn(j,i)
            bNodes(k) = 1
         enddo
      enddo

c --- check for duplicate midsides - means it is on the inside
      do i = 1,numElems
         do j = 5,8
            k = ElemConn(j,i)
            do l = i+1,numElems
               do m = 5,8
                  n=ElemConn(m,l)
                  if (k.eq.n) then
                      bNodes(k) = 0
                      goto 666
                  endif
               enddo
            enddo
 666        continue
         enddo
      enddo

c --- if either adjacent is on the border, so is the corner
      do i = 1,numElems
         k = ElemConn(1,i)
         l = ElemConn(5,i)
         m = ElemConn(8,i)
         if (bNodes(l).eq.1 .or. bNodes(m).eq.1) bNodes(k) = 1
      enddo
      do i = 1,numElems
         do j = 2,4
            k = ElemConn(j,i)
            l = ElemConn(j+3,i)
            m = ElemConn(j+4,i)
            if (bNodes(l).eq.1 .or. bNodes(m).eq.1) bNodes(k) = 1
         enddo
      enddo
c

#if defined(DEBUG)
      call SubEnd('GetBoundNodes')
#endif

      return
      end

