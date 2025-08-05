      subroutine DistortMessage (lvl)
c --- primary function: process distorted element diagnostic error messages
c
c  input arguments:
c    lvl      (int,sc,in)  - error level to be passed to erhandler (3 is for fatal)
c
#include "impcom.inc"
#include "constants.inc"

      integer lvl

#if defined(DEBUG)
      call SubBeg('DistortMessage')
#endif
      call PutErr (1)
      call erhandler('DistortMessage',5000,lvl,
     & 'Some elements have become highly distorted.',0.0d0,' ')

#if defined(DEBUG)
      call SubEnd('DistortMessage')
#endif

      return
      end
