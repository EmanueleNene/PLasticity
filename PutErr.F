      subroutine PutErr (val)
c --- primary function: update error global variabe
c
#include "impcom.inc"
#include "constants.inc"
#include "err_info.inc"
c
      integer   val

#if defined(DEBUG)
      call SubBeg('PutErr')
#endif

      call lockSMP (W_LockErr)

      E_Key = val

      call unlockSMP (W_LockErr)

      if (val.eq.1) then
          E_Key_Prev = 1
      endif
      
#if defined(DEBUG)
      call SubEnd('PutErr')
#endif
      return
      end
