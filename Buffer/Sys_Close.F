      function Sys_Close (nunit,keep_in)
c --- primary function: close a file
c
c  input arguments:
c     nunit    (int,sc,in)       - unit number
c     keep     (chr,sc,in)       - delete flag: 'K'=keep, 'D'=delete
c
c  output arguments:
c     Sys_Close   (int,sc,out)   - Sys_Close results: 0=normal, .ne.0=iostat error
c
#include "impcom.inc"
#include "constants.inc"
#include "file_sys.inc"
c
      integer  Sys_Close,nunit,i,ios
      character*(*)   keep_in
      character*1  keep
      character*8  astat(2)

      data  astat / 'KEEP    ','DELETE  ' /

#if defined(DEBUG)
      call SubBeg('Sys_Close')
#endif

c --- truncate to character*1
      keep=keep_in

      if (B_cUnit(nunit) .eq. -1) then
c --- fortran file
          i = 2
          if (keep .eq. 'K') i = 1
          ios = i
          close (unit=nunit,status=astat(i),iostat=ios)
      else
c --- c file
          call f_close  (B_cUnit(nunit),ios)
          B_cUnit(nunit) = -1
          if (keep .eq. 'D') call f_delet (ios,B_PathAndName(nunit))
          B_PathAndName(nunit) = ' '
      endif

      Sys_Close = ios

#if defined(DEBUG)
      call SubEnd('Sys_Close')
#endif

      return
      end
