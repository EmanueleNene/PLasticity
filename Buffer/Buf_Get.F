      function Buf_Get (ipage,icur,nblk)
c --- primary function: subroutine to read in a buffered I/O block
c
c  input arguments:
c      ipage   (int,sc,in)      - page number to get
c      nblk    (int,sc,in)      - buffer set number
c
c  output arguments:
c      icur    (int,sc,out)     - block number used
c
#include "impcom.inc"
#include "file_sys.inc"
#include "constants.inc"
c

      integer   Sys_Rd

      integer  Buf_Get, ipage, icur, nblk,
     &         k, i, lrec, pagenew
      integer*8  ipage8

      pointer (piBufferL,Buffer4)
      integer*4  Buffer4(*)

c --- set return value to zero
      Buf_Get = 0

c --- find available block
      icur = 0
      do i = 1,B_NumBlock(nblk)
         if (B_PageNum(i,nblk).eq.0) then
             icur = i
             goto 11
         endif
      enddo

c --- read in new block
  11  if (B_RwStat(nblk).ne.2 .or. ipage.le.B_MaxWrRec(nblk)) then

          piBufferL = B_StartPtr(nblk)
          k = B_StartWord(icur,nblk)
          lrec = B_PageLen(nblk)

          call lockSMP (W_LockBin)

          pagenew = ipage
#if defined(FLIMIT)
c --- patch for unix file length limit
          call Sys_Lim (pagenew,lrec)
#endif
c --- read in the new record
          i = Sys_Rd (B_FortUnit(nblk),pagenew,lrec,Buffer4(k+1))

          if (i.eq.-1) then
              Buf_Get = -1
              B_MaxWrRec(nblk) = ipage - 1
          elseif (i .gt. 0) then
              call FatalError(2,B_FortUnit(nblk),i)
          endif

          call unlockSMP (W_LockBin)
      endif

c --- update pointers and counters if no eof
      if (Buf_Get.ne.-1) then
          B_CurPage(nblk) = B_CurPage(nblk) + 1
          B_PageNum(icur,nblk) = ipage
          ipage8 = ipage
          B_EndLoc(icur,nblk) = ipage8 * B_PageLen(nblk)
          B_StartLoc(icur,nblk) = B_EndLoc(icur,nblk)
     &                          - B_PageLen(nblk)
          B_AgeCount = B_AgeCount + 1
          B_PageAge(icur,nblk) = B_AgeCount
          B_BlockStat(icur,nblk) = 0
      endif

      return
      end
