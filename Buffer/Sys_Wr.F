      function Sys_Wr (nunit,nrec,leng,iarray)
c --- primary function: write random record
c
c  input arguments:
c     nunit    (int,sc,in)       - file unit number
c     nrec     (int,sc,in)       - record number
c     leng     (int,sc,in)       - length of array to write (i4 words)
c     iarray   (int,ar(leng),in) - output buffer
c
c  output arguments:
c     Sys_Wr   (int,sc,out)      - 0 = success
c                                - non-zero = failure
c
#include "impcom.inc"
#include "constants.inc"
#include "file_sys.inc"
c
      integer  Sys_Wr,nunit,nrec,leng,iarray(leng),ios
      integer  lposit(2)

      lposit(1) = B_ByteInWord * leng
      lposit(2) = nrec - 1
c --- position file
      call f_lseek (B_cUnit(nunit),lposit(1),ios)
      if (ios .ne. -1) then
c --- write buffer
          call f_write (B_cUnit(nunit),B_ByteInWord*leng,iarray(1),ios)
      endif
c --- set function return code
      if (ios .eq. leng*B_ByteInWord) then
c --- successful write
          Sys_Wr = 0
      elseif (ios .lt. 0) then
c --- write error
          Sys_Wr = -ios
      else
c --- end of file error - probably full disk
          Sys_Wr = -1
      endif

      return
      end
