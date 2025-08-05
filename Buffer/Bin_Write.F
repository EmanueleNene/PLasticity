      function Bin_Write (nblk,LongLocL,nLengL,ivect4)
c --- primary function: buffer write routine
c
c  input arguments:
c     nblk      (int,sc,in)        - block number
c     LongLocL  (int*,sc,in)       - location on file
c     nLengL    (int*8,sc,in)      - record length
c     ivect4    (int*4,ar(*),in)   - vector containing record to be output
c
c  output arguments:
c     Bin_Write    (int,sc,out)    - =0: successful record write
c                                    <0: error
c
#include "impcom.inc"
#include "file_sys.inc"
c
      integer   Bin_Write, Buf_Loc

      integer*8  LongLocL, nL, jL, i1L, l1L, nLengL, kL, i
      integer    nblk, icur, i2
      integer*4  ivect4(*)


      pointer (piBufferL,Buffer4)
      integer*4  Buffer4(*)

      jL  = 0
      nL  = LongLocL
      l1L = nLengL

c --- locate proper block
  10  continue
         i1L = Buf_Loc(nL,icur,nblk)
         if (i1L .eq. -1) then
             Bin_Write = -1
             goto 999
         else
             Bin_Write = 0
         endif
         piBufferL = B_StartPtr(nblk)

         kL = l1L
         i2 = nL - B_StartLoc(icur,nblk)
         if (i2+l1L .gt. B_PageLen(nblk)) kL = B_PageLen(nblk) - i2

c --- store the vector
         do i = 1,kL
            Buffer4(i1L+i) = ivect4(jL+i)
         enddo

         B_BlockStat(icur,nblk) = 1
         if (kL .ge. l1L) goto 999

         jL = jL + kL
         l1L = l1L - kL
         nL = B_EndLoc(icur,nblk)
      goto 10

  999 continue

      return
      end
