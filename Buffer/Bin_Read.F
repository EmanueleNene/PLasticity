      function Bin_Read (nblk,LongLocL,nLengL,ivect4)
c --- primary function: buffer read routine
c
c  input arguments:
c     nblk      (int,sc,in)        - block number
c     LongLocL  (int*8,sc,in)      - current location in 32 bit word buffer
c     nLengL    (int*8,sc,in)      - record length (32 bit words)
c
c  output arguments:
c     Bin_Read  (int,sc,out)       - =0: successful record read
c                                    <0: error
c     ivect4    (int*4 ,ar(*),out) - vector containing record read
c
#include "impcom.inc"
#include "file_sys.inc"
c
      integer  Bin_Read, Buf_Loc

      integer*8  LongLocL, nL, jL, l1L, nLengL, i1L, kL, i
      integer  nblk, icur, i2
      integer*4  ivect4(*)

      pointer (piBufferL,Buffer4)
      integer*4  Buffer4(*)

      jL = 0
      nL = LongLocL
      l1L = nLengL

c --- locate proper block
  10  continue  
         i1L = Buf_Loc(nL,icur,nblk)
         if (i1L .eq. -1) then
             Bin_Read = -1
             goto 999
         else
             Bin_Read = 0
         endif
         piBufferL = B_StartPtr(nblk)

         kL = l1L
         i2 = nL - B_StartLoc(icur,nblk)
         if (i2+l1L .gt. B_PageLen(nblk)) kL = B_PageLen(nblk) - i2

c --- get the vector
         do i = 1,kL
            ivect4(jL+i) = Buffer4(i1L+i)
         enddo

         if (kL .ge. l1L) goto 999

         jL = jL + kL
         l1L = l1L - kL
         nL = B_EndLoc(icur,nblk)
      goto 10

  999 continue

      return
      end
