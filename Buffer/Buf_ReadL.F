      subroutine Buf_ReadL (nblk,LongLocL,lengL,ivect,kbfint)
c --- primary function: buffer read of a long record
c
c  input arguments:
c   nblk     (int,sc,in)      - block number
c
c  output arguments:
c   LongLocL (int*8,sc,inout) - location of the current position on the file
c   lengL    (int*8,sc,in)    - number of words to read from ivect
c   ivect    (int,ar(*),out)  - results (can be either integer or double precision)
c   kbfint   (int,sc,out)     - key for type
c                               = 0  double precision data
c                               > 0  integer data (the same as leng)
c
#include "impcom.inc"
c
      integer     nblk, kbfint, ivect(*), leng
      integer*8   LongLocL, lengL, iStopL, iStartL

      integer     RECLONG
      parameter  (RECLONG=2000000000)

      iStartL = 0
  10  iStopL = iStartL + RECLONG
          if (iStopL.gt.lengL) iStopL = lengL
          leng = iStopL - iStartL
          call Buf_Read (nblk,LongLocL,leng,ivect(iStartL+1),kbfint)
          if (leng.le.0) goto 999
          iStartL = iStopL
      if (iStopL.lt.lengL) goto 10

 999  continue
      return
      end

