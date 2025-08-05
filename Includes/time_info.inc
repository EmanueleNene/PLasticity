      integer           T_MaxLev,T_MaxSubNum
      parameter        (T_MaxLev=50,T_MaxSubNum=500)
      character*32      T_SubName
      integer           T_Stat,T_CurLev,
     &                  T_NumPerLev,T_CurSubNum
      integer*8         T_NumCalls
      double precision  T_SubTime
      common /time_info/ 
     &             T_Stat,T_CurLev, 
     &             T_SubTime(T_MaxSubNum,T_MaxLev),
     &             T_NumPerLev(T_MaxLev),
     &             T_CurSubNum(T_MaxLev),
     &             T_NumCalls(T_MaxSubNum,T_MaxLev),
     &             T_SubName(T_MaxSubNum,T_MaxLev)
