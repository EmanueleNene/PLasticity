      subroutine GetPath (extn,pname)
c --- primary function: return the full file name for a file
c
c  input arguments:
c     extn     (ch*8,sc,in)           - the standard file extension
c                                       may be lower case
c
c  output arguments:
c     pname    (ch*B_PathLen,sc,out) - the full path name of file
c
#include "impcom.inc"
#include "file_sys.inc"
c
      character*(B_ExtLen) extn
      character*(B_PathLen) pname
      integer kerr

      call MakeName (B_JobName,extn,B_Dir,pname,kerr)

      return
      end
