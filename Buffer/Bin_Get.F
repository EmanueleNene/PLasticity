      function Bin_Get (nblk)
c --- primary function: set the base offset for large file writing
c
c  input arguments:
c     nblk     (int,sc,in)         - buffer set number
c     offset   (int*8,sc,in)       - the base address for following operations
c
#include "impcom.inc"
#include "file_sys.inc"
c
      integer     nblk
      integer*8   Bin_Get

      if (nblk .gt. 0) then
          if (B_OffsetStat(nblk) .eq. 0) then
              Bin_Get = B_CurAddr(nblk)
          else
              Bin_Get = 0
          endif
      else
          Bin_Get = -BIGINT
      endif

      return
      end
