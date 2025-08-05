      subroutine DispToRST (loc8,leng,nNodes,nbuf,u,uva)
c --- primary function: write displacement vector to result file
c
c  input arguments:
c     loc8     (int*8,sc,in)     - start position of vector
c     leng     (int,sc,in)       - length of u vector (dp words)
c     nNodes   (int,sc,in)       - number of nodes in compressed form
c     nbuf     (int,sc,in)       - the open buffer for displacment reading
c     uva      (int,sc,in)       - write disp
c     u        (dp,ar(leng),in)  - the u vector
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
c
      PTRFTN   MemAlloc
      pointer (pdWorkL,dWork)
      pointer (piNodesL,Nodes)
      integer  Nodes(*)
      double precision  dWork(*)

      integer*8         loc8
      integer           leng, nNodes, nbuf, i, j, k,
     &                  uva, n, ndof
      double precision  u(leng)

#if defined(DEBUG)
      call SubBeg('DispToRST')
#endif

c --- packed form
      piNodesL = MemAlloc(nNodes,W_LenInt,'piNodesL        ')
      pdWorkL  = MemAlloc(leng,W_LenDouble,'pdWorkL         ')
      ndof = leng/nNodes
      n = 0
      j = 0
      k = 0
      do i = 1,nNodes
         if (btest(M_ElemOutMask(i),uva)) then
             k = k + 1
             Nodes(k) = i
             call v_move (u(j+1),dWork(n+1),ndof)
             n = n + ndof
         endif
         j = j + ndof
      enddo

c --- write to file
      if (n .gt. 0) then
          if (n .eq. leng) then
              j = leng*W_IntPerDp
              call Buf_Write (nbuf,loc8,j,u(1),0)
          else
              j = n*W_IntPerDp
              call Buf_Write (nbuf,loc8,j,dWork(1),0)
              call Buf_Write (nbuf,loc8,k,Nodes(1),k)
          endif
      endif

      call MemFree (pdWorkL)
      call MemFree (piNodesL)

#if defined(DEBUG)
      call SubEnd('DispToRST')
#endif

      return
      end
