      subroutine FatalError(er_type,unit,ios)
c --- primary function: to stop when a i/o error occurs
c
c  input arguments:
c     er_type  (int,sc,in)       - er_type 1: write error, 2: read error
c     unit     (int,sc,in)       - i/o unit number where error occured
c     ios      (int,sc,in)       - the "iostat" error code
c
c
#include "impcom.inc"
c
      integer    GetWriteUnit

      integer   er_type, unit, ios, iott

      iott = GetWriteUnit()

      if (er_type.eq.1) then
          write (iott,2000) unit,ios
          call StopAll(er_type)
      elseif (er_type.eq.2) then
          write (iott,2001) unit,ios
          call StopAll(er_type)
      else
          call StopAll(er_type)
      endif

      stop


 2000 format(/' *** FatalError ***'/' Write error on unit= ',
     & i6,' possible full disk.' / ' I/O status error= ',i6)
 2001 format(/' *** FatalError ***'/' Read error on unit= ',
     & i6,' possible corrupt or missing file.' /
     & ' I/O status error= ',i6)
     
      return
      end
