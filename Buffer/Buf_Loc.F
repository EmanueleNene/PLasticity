      function  Buf_Loc (nL,icur,nblk)
c --- primary function: get loc-1 of a buffer data-1
c
c  input arguments:
c     nL        (int*8,sc,in)     - location in file
c     nblk      (int,sc,in)       - block number
c
c  output arguments:
c     Buf_Loc   (int,sc,out)      - >=0, location in buffer buffer
c                                    <0, error
c     icur      (int,sc,out)      - current buffer
c
#include "impcom.inc"
#include "file_sys.inc"
c
      integer   Buf_Put, Buf_Get

      integer*8  nL, iL
      integer  curPage, icur, nblk, keof, i
      integer  Buf_Loc

      integer  nCount
      data  nCount /100/

      nCount = nCount + 1

c --- quick test for current page
      icur = B_CurBlock(nblk)
      if (icur .gt. 0) then
          iL = nL - B_StartLoc(icur,nblk)
          if (iL .ge. 0) then
              if (iL .lt. B_PageLen(nblk)) then
                  Buf_Loc = iL + B_StartWord(icur,nblk)
                  if (nCount .ge. 100) then
                      nCount = 0
                  endif
                  goto 999
              endif
          endif
      endif

c --- find page number
      curPage = ishft(nL,-B_PageSizeBit(nblk))+1   ! use shift for performance

c --- async activities
      if (nCount .ge. 100) then
          nCount = 0
      endif

c --- check for location on current page
      icur = B_CurBlock(nblk)
      if (icur .gt. 0) then
          if (B_PageNum(icur,nblk) .eq. curPage) then
              Buf_Loc = nL - B_StartLoc(icur,nblk)
     &                     + B_StartWord(icur,nblk)
            goto 999
          endif
      endif

c --- see if it is an active (in-memory) page
      icur = 0
      do i = 1,B_NumBlock(nblk)
         if (B_PageNum(i,nblk) .eq. curPage) then
             icur = i
             goto 11
         endif
      enddo

c --- swap pages
  11  keof = 0

      if (B_NumBlock(nblk) .lt. 2) then
c --- dump buffer to disc
          if (icur .lt. 1) then
              keof = Buf_Put(nblk,-1)
              keof = Buf_Get(curPage,icur,nblk)
          endif
      else
c --- multiple buffers
         if (icur .lt. 1) then
c --- data page not in memory
             if (B_CurPage(nblk) .ge. B_NumBlock(nblk)) then
                 keof = Buf_Put(nblk,curPage)
             endif
             keof = Buf_Get(curPage,icur,nblk)
         endif
      endif

c --- switch to next page
      B_CurBlock(nblk) = icur

      if (keof .ne. -1) then
          Buf_Loc = nL - B_StartLoc(icur,nblk)
     &                 + B_StartWord(icur,nblk)
      else
          Buf_Loc = -1
      endif

  999 continue
  
      return
      end
