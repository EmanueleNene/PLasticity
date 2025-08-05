      subroutine Buf_Flush (nblk)
c --- primary function: flush written records in a blocked file
c
c  input arguments:
c     nblk      (int,sc,in)       - the block number to flush
c
#include "impcom.inc"
#include "file_sys.inc"
c
      integer   Sys_Wr

      integer  nblk, i, j, ival, minval, lrec
      integer  k

      pointer (piBufferL,Buffer4)
      integer*4  Buffer4(*)


#if defined(DEBUG)
      call SubBeg('Buf_Flush')
#endif

      if (nblk .le. 0) goto 999

c --- skip if non-buffered I/O
      if (B_PageLen(nblk).eq.0) goto 999

      if (B_FortUnit(nblk) .gt. 0) then

c --- write out used pages in ascending order
          if (B_RwStat(nblk) .ne. 1) then
   5          minval = BIGINT
              ival = 0
              do i = 1,B_NumBlock(nblk)
                 if (B_BlockStat(i,nblk) .ne. 0) then
                     if (B_PageNum(i,nblk) .lt. minval) then
                         minval = B_PageNum(i,nblk)
                         ival = i
                     endif
                 endif
              enddo

              if (ival .gt. 0) then
                  lrec = B_PageLen(nblk)
                  j = B_PageNum(ival,nblk)
                  if (j .gt. B_MaxWrRec(nblk)) B_MaxWrRec(nblk) = j
                  piBufferL = B_StartPtr(nblk)
                  k = B_StartWord(ival,nblk)
                  i = Sys_Wr(B_FortUnit(nblk),j,lrec,Buffer4(k+1))
                  B_BlockStat(ival,nblk) = 0
                  goto 5
              endif
          endif
      endif

  999 continue

#if defined(DEBUG)
      call SubEnd('Buf_Flush')
#endif

      return
      end
