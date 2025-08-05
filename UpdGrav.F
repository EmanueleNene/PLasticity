      subroutine UpdGrav
c --- primary function:   compute the current value of the accelerations
c --- secondary function: set keys if any non-zero accelerations are present
c
#include "impcom.inc"
#include "solu_info.inc"
#include "acelcm.inc"
c
      integer i

c --- update accelerations at this time point
      A_Gravity = 0
      do i = 1,3
c --- refer to ANSYS common
         accel(i) = acbeg(i) + (acend(i) - acbeg(i))*A_CurLoadFact
c --- set the appropriate keys if any accelerations are present
         if (accel(i).ne.0.0d0) A_Gravity = 1
      enddo
      
      return
      end

