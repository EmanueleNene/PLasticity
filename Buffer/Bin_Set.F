      function Bin_Set (nblk,nunit,ikeyrw,npages,pname)
c --- primary function: initialize paging space for a blocked binary file.
c                       Bin_Set should be used to open a blocked file
c                       Bin_Clo should be used to close the file.
c
c  input arguments:
c     nblk      (int,sc,in)       - block number (1 to B_MaxBuf max)
c     nunit     (int,sc,in)       - fortran unit number for the file
c     ikeyrw    (int,sc,in)       - read/write flag
c                                   =  0 - both read & write
c                                   =  1 - read
c                                   =  2 - write
c       NOTE: 0 may write, but the file length may not be extended and
c               the file may or may not exist
c             1 reads only, but the file protection must set set to "rw"
c             2 may extend the file length and the file is a new file
c     npages    (int,sc,in)       - number of pages
c     pname     (chr,ar(*),in)    - name of the file
c
c  output arguments:
c     Bin_Set    (int,func,out)    -  error status
c                                  =  0 - no error
c                                  <> 0 - error occurred
c
#include "impcom.inc"
#include "constants.inc"
#include "file_sys.inc"
c
      integer  Sys_Open, PowerOf2
      PTRFTN   MemAllocL

      integer*8  nL
      integer  nblk, nunit, npages, ikeyrw, paglen, i, keyrw
      character*(*) pname
      character*4   rw
      character*16   stat(0:2), BufName
      integer  Bin_Set

c --- file existence statuses
      data  stat / 'UNKNOWN ','OLD     ','NEW     ' /


#if defined(DEBUG)
      call SubBeg('Bin_Set')
#endif
      
c      npages = min(npages,B_MaxPage)
      
      Bin_Set = 0
      paglen = PowerOf2(B_PageSize,B_PageSizeBit(nblk))

      if (nblk.le.0 .or. nblk.gt.B_MaxBuf) goto 999

c --- check for already open
      if (B_FortUnit(nblk) .gt. 0) goto 999

      rw = 'BOTH'
      keyrw = ikeyrw

      B_NumBlock(nblk) = npages
      B_FortUnit(nblk) = nunit

c --- number of 32 bit words
      B_PageLen(nblk) = paglen
      B_RwStat(nblk) = keyrw
      B_CurBlock(nblk) = 0
      B_CurPage(nblk) = 0
      B_CurAddr(nblk) = 0

      do i = 1,npages
         B_EndLoc(i,nblk) = i * B_PageLen(nblk)
         B_StartLoc(i,nblk) = B_EndLoc(i,nblk) - B_PageLen(nblk)
         B_PageAge(i,nblk) = 0
         B_BlockStat(i,nblk) = 0
         B_StartWord(i,nblk) = B_StartLoc(i,nblk)
         B_PageNum(i,nblk) = 0
      enddo

c --- get memory for file buffer
      if (paglen.gt.0) then
          write (BufName,1001) nblk
 1001     format ('Buffer',i2)
          nL = paglen
          nL = nL*npages
#if defined (FULL_64)
          nL = nL/2
#endif
          B_StartPtr(nblk) = MemAllocL(nL,W_LenInt,BufName)
      else
          call MemFree (B_StartPtr(nblk))
      endif

      if (keyrw .ne. 2) then
          do i = 1,npages
             B_PageNum(i,nblk) = 0
          enddo
      endif

c --- open random access file
      if (nunit .gt. 0) then
         Bin_Set = Sys_Open (nunit,pname,stat(keyrw),'DIRECT',
     &                       'UNFORM',rw,'BEGIN')

         if (Bin_Set .ne. 0) then
             B_FortUnit(nblk) = 0
             call MemFree (B_StartPtr(nblk))
             goto 999
         endif

      endif

      B_OpenStat(nblk) = 1
      B_CurPos(nblk) = 0
      B_MaxLen(nblk) = 0
      B_MaxWrRec(nblk) = 0


  999 continue

#if defined(DEBUG)
      call SubEnd('Bin_Set')
#endif

      return
      end
