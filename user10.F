      function user10()
c --- primary function: reset random seed
c
#include "impcom.inc"
#include "ansysdef.inc"
#include "maccom.inc"
c
      integer user10
      double precision v(1)

      user10 = 0

      mcseed = 2.0d0
      call v_rand (v(1),1,mcseed) 
      
      return
      end