      function Buf_Inq (nblk)
c --- primary function: get data about a block i/o buffer status
c
c  input arguments:
c     nblk      (int,sc,in)        - the block number for the inquiry
c
c  output arguments:
c     Buf_Inq    (int,func,out)    - the returned status
c
#include "impcom.inc"
#include "file_sys.inc"
c
      integer*8  Buf_Inq, jL
      integer  nblk

      jL = 0

      if ((nblk.gt.0) .and. (nblk.le.B_MaxBuf)) then
          jL = B_RwStat(nblk)
      endif

      Buf_Inq = jL

      return
      end
