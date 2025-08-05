      subroutine GeomToRST (eorder)
c --- primary function: put geometry into the result file
c
c  input arguments:
c    eorder    (int,ar(*),in)    - element order vector
c
#include "impcom.inc"
#include "constants.inc"
#include "RST_file.inc"
#include "dof_info.inc"
c
      integer   etyget, etyiqr, ndgall, elmget, rlinqr, rlget
      PTRFTN    MemAlloc, MemAllocL

      integer*8  loc8, leng8, idxlocL, jdxlocL, iloc8, iOffsetL
      integer  eorder(D_NumEl), ivect(40), n, kbf, maxr,
     &         i, j , maxet, k, ityp, elChar(IELCSZ),
     &         nodes(W_NdMax+EL_DIM),LengGeomHead


      pointer (plWorkL,lWork)
      pointer (piWorkL,iWork)
      pointer (pdWorkL,dWork)      
      pointer (pdNodeL,dNode)
      integer*8  lWork(*)
      integer    iWork(*)
      double precision  dWork(*), dNode(7,*)


#if defined(DEBUG)
      call SubBeg('GeomToRST')
#endif

c --- read header record
      iloc8 = W_HedLen
      loc8 = 0
      call Bin_Put (R_BlockNum,loc8)
      n = 40
      call Buf_Read (R_BlockNum,iloc8,n,ivect(1),kbf)

c --- store geometry index location
      ivect(16) = R_FilePos
      iloc8 = W_HedLen
      n = 40
      call Buf_Write (R_BlockNum,iloc8,n,ivect(1),n)

c --- set up the geometry index
      LengGeomHead = 40
      call vi_zero (ivect(1),LengGeomHead)
      idxlocL = R_FilePos
      R_FilePos = R_FilePos + LengGeomHead + 3

c --- store element types
      maxet = etyiqr(0,DB_MAXDEFINED)
      ivect(2) = maxet
      ivect(19) = IELCSZ
      if (maxet .gt. 0) then
          jdxlocL  = R_FilePos
          call MakeShortInt (R_FilePos,ivect(21),ivect(22))
          iOffsetL = jdxlocL
          piWorkL = MemAlloc(maxet,W_LenInt,'piWorkL        ')
          call vi_zero (iWork(1),maxet)
          R_FilePos = R_FilePos + maxet + 3  ! space for relative index
          do i = 1,maxet
             n = etyget(i,elChar(1))
             if (n .ne. 0) then
                 n = abs(n)
                 iWork(i) = R_FilePos - iOffsetL
                 call Buf_Write (R_BlockNum,R_FilePos,n,elChar(1),n)
             endif
          enddo
          call Buf_Write (R_BlockNum,jdxlocL,maxet,iWork(1),maxet)
          call MemFree (piWorkL)
      endif

c --- store the real constants (if those are read in, we do not need RESUME in post-processor)
      maxr = rlinqr(0,DB_MAXDEFINED)
      ivect(3) = maxr
      ivect(20) = rlinqr(0,DB_MAXRECLENG)
      if (maxr .gt. 0) then
          jdxlocL = R_FilePos
          call MakeShortInt (R_FilePos,ivect(23),ivect(24))
          iOffsetL = jdxlocL
          n = rlinqr(0,DB_MAXRECLENG)
          piWorkL = MemAlloc(maxr,W_LenInt,'piWorkL        ')
          pdWorkL = MemAlloc(n,W_LenDouble,'pdWorkL        ')
          call vi_zero (iWork(1),maxr)
          R_FilePos = R_FilePos + maxr + 3  ! space for relative index
          do i = 1,maxr
             n = rlget(i,dWork(1))
             if (n .ne. 0) then
                 n = abs(n)
                 iWork(i) = R_FilePos - iOffsetL
c --- save real constants
                 call Buf_Write (R_BlockNum,R_FilePos,n*W_IntPerDp,
     &                           dWork(1),0)
             endif
          enddo
c --- save the index of the real constants in the file
          call Buf_Write (R_BlockNum,jdxlocL,maxr,iWork(1),maxr)
          call MemFree (pdWorkL)
          call MemFree (piWorkL)
      endif

c --- store the nodes
      n = D_NumNd
      ivect(4) = n
      if (n .gt. 0) then
          call MakeShortInt (R_FilePos,ivect(27),ivect(28))
          iloc8 = n
          iloc8 = 7*iloc8 + 7
          pdNodeL = MemAllocL(iloc8,W_LenDouble,'NodeList       ')
           do j = 1,n
              do i = 1,7
                 dNode(i,j) = 0.0d0
              enddo
           enddo
           j = 1
           do i = 1,n
              k = ndgall(i,dNode(2,j))
              if (k .ne. 0) then
                  dNode(1,j) = i
                  j = j + 1
              endif
           enddo
           do i = 1,n
              call Buf_Write (R_BlockNum,R_FilePos,
     &                        7*W_IntPerDp,dNode(1,i),0)
           enddo
           call MemFree (pdNodeL)
      endif

c --- store the elements
      ivect(5) = D_NumEl
      ivect(18) = W_NdMax
      jdxlocL = R_FilePos
      call MakeShortInt (R_FilePos,ivect(29),ivect(30))
      iOffsetL = jdxlocL
      plWorkL = MemAlloc(D_NumEl,W_LenLongInt,'plWorkL        ')
      do i = 1,D_NumEl
         lWork(i) = 0
      enddo
      R_FilePos = R_FilePos + W_IntPerLong*D_NumEl + 3  ! storage for index
      do i = 1,D_NumEl
         k = eorder(i)
         if (k .gt. 0) then
             call vi_zero (nodes(EL_DIM+1),W_NdMax)
             n = elmget(k,nodes(1),nodes(EL_DIM+1))
c --- store element number in 4-th position
c     note, if date is stored, it has to be 9-th position
             nodes(4) = k
             if (n .ne. 0) then
                 ityp = nodes(EL_TYPE)
                 n = EL_DIM + abs(n)
                 lWork(i) = R_FilePos - iOffsetL
                 call Buf_Write (R_BlockNum,R_FilePos,n,nodes(1),n)
             endif
         endif
      enddo

      i = W_IntPerLong*D_NumEl
      call MakeLongIntV (D_NumEl,lWork(1))
      call Buf_Write (R_BlockNum,jdxlocL,i,lWork(1),i)
      call MemFree (plWorkL)

c --- write the updated index
      call Buf_Write (R_BlockNum,idxlocL,LengGeomHead,
     &                ivect(1),LengGeomHead)

c --- update total file length
      iloc8 = W_HedLen
      n = 40
      call Buf_Read (R_BlockNum,iloc8,n,ivect(1),kbf)
      ivect(10) = R_FilePos
      leng8 = R_FilePos
      call MakeShortInt (leng8,ivect(23),ivect(24))
      iloc8 = W_HedLen
      n  = 40
      call Buf_Write (R_BlockNum,iloc8,n,ivect(1),n)


#if defined(DEBUG)
      call SubEnd('GeomToRST')
#endif

      return
      end
