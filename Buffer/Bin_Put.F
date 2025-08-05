      subroutine Bin_Put (nblk,offset)
c --- primary function: set the base offset for large file writing
c
c  input arguments:
c     nblk     (int,sc,in)         - buffer set number
c     offset   (integer*8,sc,in)   - the base address for following operations
c
#include "impcom.inc"
#include "file_sys.inc"
c
      integer     nblk
      integer*8   offset

#if defined(DEBUG)
      call SubBeg('Bin_Put')
#endif

      if (nblk .gt. 0) then
          if (B_OffsetStat(nblk) .eq. 0) then
              B_CurAddr(nblk) = offset
          else
              B_CurAddr(nblk) = 0
              B_OffsetStat(nblk) = offset
              if (offset .le. 0) B_OffsetStat(nblk) = 1
          endif
      endif


#if defined(DEBUG)
      call SubEnd('Bin_Put')
#endif

      return
      end
