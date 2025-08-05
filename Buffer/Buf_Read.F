      subroutine Buf_Read (nblk,LongLocL,leng,ivect,kbfint)
c --- primary function: buffer read
c
c  input arguments:
c   nblk     (int,sc,in)     - block number
c
c  output arguments:
c   LongLocL (int*8,sc,inout)- location of the current position on the file
c   leng     (int,sc,inout)  - number of items read zero in case of a problem
c   ivect    (int,ar(*),out) - results (can be either integer or double precision)
c   kbfint   (int,sc,out)    - key for type
c                              = 0  double precision data
c                              > 0  integer data (the same as leng)
c
#include "impcom.inc"
#include "file_sys.inc"
c
      integer   Bin_Read
      integer*8 MakeLongInt
      integer*8 LongLocL, loc8, nL
      integer   nblk, leng, ivect(*), kbfint, i,
     &          iwork(2), lengin
      integer*4 iwork4(2)

      if (nblk .le. 0) goto 999

      if (LongLocL .eq. -1) then
          loc8 = B_CurPos(nblk)
      else
          loc8 = LongLocL + B_CurAddr(nblk)
      endif

      lengin = leng
      leng = 0
      if (B_FortUnit(nblk) .le. 0) goto 999

c --- read initial two 32 bit words
      nL = 2
      i = Bin_Read(nblk,loc8,nL,iwork4(1))
      if (i .eq. -1) goto 900

      if (B_PageFormat(nblk).eq.0) then
          iwork(1) = iwork4(1)
          iwork(2) = iwork4(2)
      else
          call Ex_InInt (iwork4(1),iwork(1),2)
      endif

c --- large format
      if (btest(iwork(2),31)) then
          kbfint = 1
          iwork(2) = ibclr(iwork(2),31)
      else
          kbfint = 0
      endif
      if (iwork(2) .ne. 0) then
          leng = MakeLongInt(iwork(1),iwork(2))
      else
          leng = iwork(1)
      endif
      if (kbfint .gt. 0) kbfint = leng

      loc8 = loc8 + 2

c --- end-of-file condition (leng = -1)
      if (leng .le. 0) then
          loc8 = loc8 + 1
          goto 900
      endif
c --- read the data
      if (lengin .gt. 0) then
          nL = leng
          i = Bin_Read(nblk,loc8,nL,ivect(1))

          if (i .eq. -1) goto 900
          if (B_PageFormat(nblk).ne.0) then
              if (kbfint .eq. 0) then
                  call Ex_InDP (ivect(1),ivect(1),leng/2)
              else
                  call Ex_InInt (ivect(1),ivect(1),leng)
              endif
          endif
 
      endif
      loc8 = loc8 + leng

c --- skip the trailing word
      loc8 = loc8 + 1
      B_CurPos(nblk) = loc8

      goto 990

c --- end of file
  900 leng = 0

  990 LongLocL = loc8 - B_CurAddr(nblk)

  999 continue

      return
      end
