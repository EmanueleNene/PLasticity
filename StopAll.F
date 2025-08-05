      subroutine StopAll (icode)
c --- primary function: stop run
c
c  input arguments:
c     icode    (int,sc,in)       - stop error code
c                                  0 - normal exit
c
#include "impcom.inc"
c
      integer   icode, i, FormOutput

      if (icode.ge.0) then
          i = FormOutput()
      else
          stop
      endif

c --- some logic should be added here of icode.ne.0
            
      return
      end
