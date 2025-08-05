      function Buf_Put (nblk,pagehold)
c --- primary function: subroutine to write out oldest I/O block
c
c  input arguments:
c      nblk     (int,sc,in)      - buffer set number
c      pagehold (int,sc,in)      - current page number
c
c  output arguments:
c      Buf_Put  (int,sc,out)     - 0, successful
c
#include "impcom.inc"
#include "constants.inc"
#include "file_sys.inc"
c

      integer   Sys_Wr

      integer  Buf_Put, icur, nblk, pagehold,
     &         k, i, j, lrec
      double precision  dperr(2)


      pointer (piBufferL,Buffer4)
      integer*4  Buffer4(*)


c --- set return value to zero
      Buf_Put = 0
      icur = 0
      if (B_CurPage(nblk).le.0) goto 999

c --- find oldest block
      do i = 1,B_NumBlock(nblk)
         j = B_PageNum(i,nblk)
         if (j .ne. pagehold) then
             if (icur .le. 0) then
                 icur = i
             else
                 if (B_PageAge(i,nblk).lt.B_PageAge(icur,nblk)) icur = i
            endif
         endif
      enddo

      if (icur.le.0) goto 999

c --- output oldest block
      if (B_BlockStat(icur,nblk) .gt. 0) then

          j = B_PageNum(icur,nblk)
          k = B_StartWord(icur,nblk)
          piBufferL = B_StartPtr(nblk)
          call lockSMP (W_LockBin)
          lrec = B_PageLen(nblk)

c --- save maximum record number
          if (j.gt.B_MaxWrRec(nblk)) B_MaxWrRec(nblk) = j
#if defined(FLIMIT)
c --- patch for unix file length limit
          call Sys_Lim(j,lrec)
#endif
c --- write out the record
          i = Sys_Wr (B_FortUnit(nblk),j,lrec,Buffer4(k+1))

c --- write out error diagnostics
          if (i.ne.0) then
              dperr(1) = dble (B_FortUnit(nblk))
              dperr(2) = dble (i)
              call erhandler('Buf_Put',5000,2,
     &                       'I/O error. File = %i
     &                       %/The system I/O error = %i.',
     &                       dperr(1),' ')

              call FatalError(1,B_FortUnit(nblk),i)
          endif
          call unlockSMP (W_LockBin)
      endif

      B_CurPage(nblk) = B_CurPage(nblk) - 1
      B_PageNum(icur,nblk) = 0
      B_EndLoc(icur,nblk) = 0
      B_StartLoc(icur,nblk) = 0
      B_PageAge(icur,nblk) = 0
      B_BlockStat(icur,nblk) = 0

  999 continue

      return
      end
