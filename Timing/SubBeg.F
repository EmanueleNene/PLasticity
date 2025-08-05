      subroutine SubBeg (subname)
c --- primary function: mark beginning of a subroutine
c
c  input arguments:
c     subname   (char,sc,in)    - name of subroutine being entered
c
#include "impcom.inc"
      character*(*)     subname
      
#if defined(ATRACK)
      call TrackBegin (subname)
#else

#include "constants.inc"
#include "time_info.inc"
#include "solu_info.inc"
c
      integer           iProcSMP, GetWriteUnit
      integer           i, n, lev, iott
      double precision  cpu,time

      call KCL_Time(cpu,time)

      if (A_Debug .eq. -1) then
          iott = GetWriteUnit()
          call lockSMP (W_LockText)
          write (iott,*) 'Enter "',subname,'" at proc num =', iProcSMP()
          call unlockSMP (W_LockText)
      endif

c --- track only master thread now
      if (iProcSMP().gt.0 .or. T_Stat.eq.0) goto 999
      
      if (T_CurLev.gt.0) then
          n = T_CurSubNum(T_CurLev)
c --- suspend previous level
          T_SubTime(n,T_CurLev) = T_SubTime(n,T_CurLev) + cpu
      endif
      
      lev = T_CurLev + 1
c --- check for level overflow
      if (lev .gt. T_MaxLev) then
          iott = GetWriteUnit()
          call lockSMP (W_LockText)
          write (iott,1000)
 1000     format (/'Counter exceeded max level. Should never happen!')
          call unlockSMP (W_LockText)
          goto 999
      endif
c --- increment level counter
      T_CurLev = lev

c --- check previous name first
      n = T_CurSubNum(lev)
      if (n .gt. 0) then
         if (subname .eq. T_SubName(n,lev)) goto 20
      endif

c --- find another at this level
      n = T_NumPerLev(lev)
      do i = n,1,-1         ! search in reverse order
         if (subname .eq. T_SubName(i,lev)) then
             n = i
             goto 20
         endif
      enddo

c --- add new name at this level
      if (n .lt. T_MaxSubNum) n = n + 1
      T_SubName(n,lev) = subname
      T_NumPerLev(lev) = n

  20  T_NumCalls(n,lev) = T_NumCalls(n,lev) + 1

      T_CurSubNum(lev) = n
      call KCL_Time(cpu,time)
      T_SubTime(n,lev) = T_SubTime(n,lev) - cpu

 999  continue
#endif 
      return
      end
