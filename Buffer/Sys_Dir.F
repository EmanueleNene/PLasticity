      function Sys_Dir (fdir)
c --- primary function: get the path name for the user directory
c
c  output arguments:
c     fdir      (chr,sc,out)     - character string for the directory name
c     Sys_Dir   (int,sc,out)     - 0 = success
c                                - 1 = failure
c
#if defined (PCWINNT_SYS)
#if defined (__IFL) || defined (__EFL)
#if (__INTEL_COMPILER >= 800)
      USE IFPORT
#else
      USE IFLPORT
#endif
#else
      USE DFLIB
#endif
#endif
#include "impcom.inc"
#include "file_sys.inc"

      integer  Sys_Dir,ij
      character*(B_DirLen)  fdir

#if defined (PCWINNT_SYS)
      integer  length
      character*(B_DirLen)  filnm1
      character*2  fdrive

      fdir = '                                '
      fdrive = '  '
      Sys_Dir = 0
      length = getenvqq('HOMEPATH',fdir)
      length = getenvqq('HOMEDRIVE',fdrive)
      if (fdir(2:2) .ne. ' ') then
          if (fdrive .ne. '  ') then
              filnm1 = fdrive // '\' // fdir
              fdir = filnm1
          endif
          Sys_Dir=0
      else
          if (fdir(1:1) .eq. '\') then
              filnm1 = fdrive // fdir
              fdir = filnm1
          else
              Sys_Dir = 1
          endif
      endif

#else
      fdir = '                                '
      Sys_Dir = 0
      call getenv('HOME',fdir)
      if (fdir(2:2) .ne. ' ') then
          do 9000 ij = 3,B_DirLen
             if (fdir(ij:ij) .eq. ' ') then
                fdir(ij:ij) = '/'
                Sys_Dir = 0
                goto 999
             endif
 9000     continue
          Sys_Dir = 1
      endif
#endif
  999 continue

      return
      end
