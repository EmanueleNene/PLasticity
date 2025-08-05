      subroutine Buf_Write (nblk,LongLocL,leng,ivect,kbfint)
c --- primary function: buffer write
c
c  input arguments:
c   nblk   (int,sc,in)        - block number
c   leng   (int,sc,in)        - number of words to read from ivect
c   kbfint (int,sc,in)        - key for type
c                               = 0  double precision data
c                               > 0  integer data (same as leng)
c
c  output arguments:
c   LongLocL (int*8,sc,inout) - location in integer words of the current
c                               position on the file
c   ivect   (int,ar(*),out)   - vector containing record to be written
c
#include "impcom.inc"
#include "file_sys.inc"
c
      integer    SZE
      parameter (SZE=1024)
      integer*4  iwork4(SZE)

      integer    Bin_Write

      integer*8  LongLocL, loc8, iL
      integer    nblk, leng, ivect(*), kbfint, i, n, iwork(2)

      if (nblk .le. 0) goto 999

      loc8 = LongLocL + B_CurAddr(nblk)

c --- determine the number of integer*4 words to put on file
      if (kbfint .gt. 0) then
         n = leng
      else
#if defined (FULL_64)
         n = leng*2
#else
         n = leng
#endif
      endif

      if (B_FortUnit(nblk) .gt. 0) then

          if (LongLocL .eq. -1) loc8 = B_CurPos(nblk)

c --- write initial two words

          iwork(1) = n
          iwork(2) = 0

          if (kbfint .gt. 0) iwork(2) = ibset(iwork(2),31)

          if (B_PageFormat(nblk).eq.0) then
              iwork4(1) = iwork(1)
              iwork4(2) = iwork(2)
          else
              call In_ExInt (iwork(1),iwork4(1),2)
          endif

          iL = 2
          i =  Bin_Write(nblk,loc8,iL,iwork4(1))

          loc8 = loc8 + 2

c --- write the data
          iL = n
          if (B_PageFormat(nblk).eq.0) then
              i = Bin_Write(nblk,loc8,iL,ivect(1))
          else
              if (leng .gt. SZE) then
                  if (kbfint .eq. 0) then
                      call In_ExDP (ivect(1),ivect(1),leng/2)
                  else
                      call In_ExInt (ivect(1),ivect(1),leng)
                  endif
                  i = Bin_Write(nblk,loc8,iL,ivect(1))
                  if (kbfint .eq. 0) then
                      call Ex_InDP (ivect(1),ivect(1),leng/2)
                  else
                      call Ex_InInt (ivect(1),ivect(1),leng)
                  endif
              else
c --- for smaller chunks, use a temporary buffer
                  if (kbfint .eq. 0) then
                      call In_ExDP (ivect(1),iwork4(1),leng/2)
                  else
                      call In_ExInt (ivect(1),iwork4(1),leng)
                  endif
                  i = Bin_Write(nblk,loc8,iL,iwork4(1))
              endif
          endif
          loc8 = loc8 + n

c --- write the trailing word
          iL = 1
          i = Bin_Write(nblk,loc8,iL,iwork4(1))
          loc8 = loc8 + 1

      else
c --- bit bucket
          loc8 = loc8 + n + 3
      endif

      B_CurPos(nblk) = loc8
      B_MaxLen(nblk) = max(loc8,B_MaxLen(nblk))
      LongLocL = loc8 - B_CurAddr(nblk)

  999 continue

      return
      end
