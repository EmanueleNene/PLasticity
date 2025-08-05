      subroutine Buf_WriteL (nblk,LongLocL,lengL,ivect,kbfint)
c --- primary function: buffer write of a long record
c
c  input arguments:
c   nblk   (int,sc,in)        - block number
c   lengL  (int*8,sc,in)      - number of words to read from ivect
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
c
      integer     nblk, kbfint, ivect(*), leng
      integer*8   LongLocL, lengL, iStartL, iStopL

      integer     RECLONG
      parameter  (RECLONG=2000000000)

      iStartL = 0

  10  iStopL = iStartL + RECLONG
          if (iStopL.gt.lengL) iStopL = lengL
          leng = iStopL - iStartL
          call Buf_Write (nblk,LongLocL,leng,ivect(iStartL+1),kbfint)
          if (leng.le.0) goto 999
          iStartL = iStopL
      if (iStopL.lt.lengL) goto 10

 999  continue
      return
      end
         