      subroutine Buf_Inf (nblk,ival)
c --- primary function: define information for buffers status
c
c  input arguments:
c     nblk     (int,sc,in)       - the block number
c     ival     (int,sc,in)       - the new value
c
#include "impcom.inc"
#include "file_sys.inc"
c
      integer  nblk,ival

      if (nblk .gt. 0) then
          B_RwStat(nblk) = ival
      endif

      return
      end
