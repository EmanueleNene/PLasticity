      subroutine TimeIni
c --- primary function: initialize timing
c
#include "impcom.inc"
#include "time_info.inc"
c
      integer i,j
      
      T_CurLev = 0
      T_Stat = 0

      call v_zero (T_SubTime(1,1),T_MaxLev*T_MaxSubNum)
      
      do j = 1,T_MaxLev
         do i = 1,T_MaxSubNum
            T_NumCalls(i,j) = 0
         enddo
      enddo      
      
      call vi_zero(T_NumPerLev(1),T_MaxLev)
      
      T_Stat = 1   ! start timing
      
      return
      end
