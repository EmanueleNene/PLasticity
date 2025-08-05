      function GetWriteUnit ()
c --- primary function: obtain information about output
c
c  output arguments:
c     GetWriteUnit   (int,func,out)    - the output unit number
c
#include "impcom.inc"
c
      integer  GetWriteUnit, wrinqr

      GetWriteUnit = wrinqr (2)

      return
      end
