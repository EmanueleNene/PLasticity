      subroutine SubEnd (subname)
c --- primary function: mark end of a subroutine
c
c  input arguments:
c     subname     (char,sc,in)  - name of subroutine being left
c
#include "impcom.inc"
      character*(*)     subname
      
#if defined(ATRACK)
      call TrackEnd (subname)
#else

#include "constants.inc"
#include "time_info.inc"
#include "solu_info.inc"
c
      integer           iProcSMP, GetWriteUnit
      integer           n, lev, iott
      double precision  cpu, time

      call KCL_Time(cpu,time)

      if (A_Debug .eq. -1) then
          iott = GetWriteUnit()
          call lockSMP (W_LockText)
          write (iott,*) 'Leave "',subname,'" at proc num =', iProcSMP()
          call unlockSMP (W_LockText)
      endif   
      
c --- track only master thread now
      if (iProcSMP().gt.0 .or. T_Stat.eq.0) goto 999

c --- decrement level counter
      lev = T_CurLev
      T_CurLev = T_CurLev - 1

      if (T_CurLev .lt. 0) then
          iott = GetWriteUnit()
          call lockSMP (W_LockText)
          write (iott,1000)
 1000     format (/'Level counter went negative. Should never happen!')
          call unlockSMP (W_LockText)
          T_CurLev = 0
          lev = 1
      endif

c --- update timing tables
      n = T_CurSubNum(lev)
      if (n .gt. 0) T_SubTime(n,lev) = T_SubTime(n,lev) + cpu
      
      if (T_CurLev.gt.0) then
          n = T_CurSubNum(T_CurLev)
c --- resume previous level
          if (n .gt. 0) then
              call KCL_Time(cpu,time)
              T_SubTime(n,T_CurLev) = T_SubTime(n,T_CurLev) - cpu
          endif
      endif
 
 999  continue
 
#endif 
      return
      end
