      function GetErr()
c --- primary function: obtain information from errors common
c
c  output arguments:
c     GetErr   (int,sc,out)      - value corresponding to key
c
#include "impcom.inc"
#include "err_info.inc"
c
      integer   GetErr

      GetErr = E_Key

      return
      end
