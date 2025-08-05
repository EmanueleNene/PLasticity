      subroutine Buf_Ini
c --- primary function: initialize buffered system
c
#include "impcom.inc"
#include "file_sys.inc"
c
      integer  i, n, PowerOf2

#if defined(DEBUG)
      call SubBeg('Buf_Ini')
#endif

c --- initialize variables
c      B_NumOfPages = B_MaxPage
c      B_UsedPageLen = 65436

      do i = 1,B_MaxBuf
         B_FortUnit(i) = 0
      enddo
      B_AgeCount = 0
      B_PageSize = PowerOf2(B_UsedPageLen,n)
      do i = 1,B_MaxBuf
         B_MaxLen(i) = 0
         B_CurAddr(i) = 0
         B_OffsetStat(i) = 0
         B_MaxWrRec(i) = 0
         B_PageSizeBit(i) = n
         B_StartPtr(i) = PTRFTNNULL
         B_PageFormat(i) = 0
      enddo

      do i = 1,B_MaxUnit
c --- initialize i/o descriptors
         B_cUnit(i) = -1
c --- initialize path names
         B_PathAndName(i) = ' '
      enddo

#if defined(DEBUG)
      call SubEnd('Buf_Ini')
#endif

      return
      end
