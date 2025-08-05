      subroutine Bin_Clo (nblk,pstat)
c --- primary function: close a blocked file
c
c  input arguments:
c     nblk      (int,sc,in)       - the block number to close
c
c     pstat     (chr,sc,in)       - keep or delete flag
c                                   = 'K' - keep file
c                                   = 'D' - delete file
c
#include "impcom.inc"
#include "file_sys.inc"
c
      integer   Bin_Write, Sys_Close

      integer*8  bioLeng8, loc8, nLengL, iL
      integer    nblk, i
      integer*4  iwork4(1)
      integer    iwork(1)
      character  pstat
      character*1  stat, kstat

#if defined(DEBUG)
      call SubBeg('Bin_Clo')
#endif

      if (nblk.le.0 .or. nblk.gt.B_MaxBuf) goto 999
      stat = pstat

      if (B_FortUnit(nblk) .gt. 0) then
          kstat = stat
c --- check for write or read/write
          if (B_RwStat(nblk).ne.1) then
c --- check to see if written file is empty - if empty, delete
              if (B_MaxLen(nblk).eq.0) then
                  kstat = 'D'
              else
c --- rewrite header with new lengths
                  if (stat.eq.'K') 
     &                call Bin_Head(nblk,B_FortUnit(nblk),iL)
              endif
              bioLeng8 = B_PageLen(nblk)
              if (kstat.eq.'K' .and. bioLeng8.gt.0) then
c --- put -1 in last word of partial record
                  loc8 = B_MaxLen(nblk)
                  if (mod(loc8,bioLeng8).ne.0) then
                      if (B_PageFormat(nblk).eq.0) then
                          iwork4(1) = -1
                      else
                          iwork(1) = -1
                          call In_ExInt(iwork(1),iwork4(1),1)
                      endif
                      nLengL = 1
                      i = Bin_Write(nblk,loc8,nLengL,iwork4(1))
                  endif
c --- write out used pages
                  call Buf_Flush (nblk)
              endif
          endif

c --- close the file
          i = Sys_Close (B_FortUnit(nblk),kstat)

          B_FortUnit(nblk) = 0
          B_CurAddr(nblk) = 0
          B_OffsetStat(nblk) = 0
      endif

c --- release the buffer space
      call MemFree (B_StartPtr(nblk))

      B_OpenStat(nblk) = 0
      nblk = 0

  999 continue

#if defined(DEBUG)
      call SubEnd('Bin_Clo')
#endif

      return
      end
