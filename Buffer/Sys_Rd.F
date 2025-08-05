      function Sys_Rd (nunit,nrec,leng,iarray)
c --- primary function: read random record
c
c  input arguments:
c     nunit    (int,sc,in)       - file unit number
c     nrec     (int,sc,in)       - record number
c     leng     (int,sc,in)       - length of buffer to read (i4 words)
c
c  output arguments:
c     iarray   (int,ar(leng),out)- buffer being read
c     Sys_Rd   (int,sc,out)      - 0 = success
c                                - 1 = end-of-file
c                                  other than 0 or -1 = failure
c
#include "impcom.inc"
#include "constants.inc"
#include "file_sys.inc"
c
      integer  Sys_Rd,nunit,nrec,leng,iarray(leng),ios
      integer  lposit(2)

      lposit(1) = B_ByteInWord * leng
      lposit(2) = nrec - 1
c --- position file
      call f_lseek (B_cUnit(nunit),lposit(1),ios)
      if (ios .ne. -1)  then
c --- read file
          call f_read (B_cUnit(nunit),B_ByteInWord*leng,iarray(1),ios)
      endif
c --- set function return code
      if (ios .eq. 0) then
c --- end of file
          Sys_Rd = -1
      elseif (ios .lt. 0) then
c --- read error
          Sys_Rd = ios
      elseif (ios .gt. 0) then
c --- successful read
          Sys_Rd = 0
      endif

      return
      end
