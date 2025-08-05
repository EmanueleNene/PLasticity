      function Sys_Inq (fname,inqr_in)
c --- primary function: do a file system inquire
c
c  input arguments:
c    nunit    (int,sc,in)       - fortran unit number (used only for inqr='O')
c    fname    (chr,sc,in)       - name of file
c    inqr_in  (chr,sc,in)       - character key for information requested
c                                 = 'E' - return whether file exists
c                                     Sys_Inq = 1 - file exists
c                                             = 0 - file does not exist
c                                             < 0 - error occured
c                                 = 'O' - return whether file is open
c                                     Sys_Inq = 1 - file is open
c                                             = 0 - file is closed
c                                             < 0 - error occured
c                                 = 'N' - return unit number of file
c                                     Sys_Inq > 0 - unit number for file
c                                             = 0 - file not assigned to a unit
c                                             < 0 - error occured
c
c  output arguments:
c    Sys_Inq   (int,func,out)    - the returned value
c
#include "impcom.inc"
#include "constants.inc"
#include "file_sys.inc"
c
      character*(*) inqr_in
      integer  Sys_Inq, lname, ios, val
      integer  i, GetStringLen, lll
      character*1  inqr
      character*(*)  fname
      logical  lval

#if defined(DEBUG)
      call SubBeg('Sys_Inq')
#endif

      Sys_Inq = 0
      ios = 0
      lname = GetStringLen(fname)
      inqr = inqr_in(1:1)
      lval = .false.
c --- get length of name
      lll = GetStringLen (fname)
      if (inqr .eq. 'E') then
          inquire (file=fname(1:lll  ),err=90,iostat=ios,exist=lval)
          if (lval) Sys_Inq = 1
      elseif (inqr .eq. 'O') then
c --- check for c file open
          do 10 i = 1, B_MaxUnit
             if (fname(1:lll) .eq. B_PathAndName(i)(1:lll)) then
                 Sys_Inq = 1
                 goto 999
             endif
   10     continue
c --- inquire statements with open option react unpredictably
c --- for files that do not exist
          inquire (file=fname(1:lll),err=90,iostat=ios,exist=lval)
          if (lval) then
              inquire(file=fname(1:lll),err=90,iostat=ios,opened=lval)
              if (lval) Sys_Inq = 1
          else
              Sys_Inq = 0
          endif
      elseif (inqr .eq. 'N') then
c --- get unit number for attached file name
c --- check for c file
          do 20 i = 1, B_MaxUnit
             if (fname(1:lll) .eq. B_PathAndName(i)(1:lll)) then
                 Sys_Inq = i
                 goto 999
             endif
   20     continue
c --- else check for fortran open
          val = 0
          inquire(file=fname(1:lname),err=90,iostat=ios,number=val)
          Sys_Inq = val
      endif
      goto 999
  90  Sys_Inq = - abs(ios)
  999 continue

#if defined(DEBUG)
      call SubEnd('Sys_Inq')
#endif

      return
      end
