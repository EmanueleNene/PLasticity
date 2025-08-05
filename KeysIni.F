      subroutine KeysIni
c --- primary function: initialize the element formation keys
c
#include "impcom.inc"
#include "solu_info.inc"
c
#if defined(DEBUG)
      call SubBeg('KeysIni')
#endif

c --- initialize forming keys:
      call vi_zero (A_MatrReq(1),10)

      A_KeyMatUpd = 0
      A_KeyBCsUpd  = 0
      A_KeyMoistUpd = 0
      A_KeyFirstIter = 0
      A_KeyOut = 0

c --- base saved variables are always written!
      A_KeySvrUpd = 1

#if defined(DEBUG)
      call SubEnd('KeysIni')
#endif

      return
      end

