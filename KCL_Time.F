      subroutine KCL_Time (cpu,time)
c --- primary function: very precise timing
c
c  output arguments:
c     cpu       (dp,sc,out)   - cpu seconds
c     time      (dp,sc,out)   - wall seconds
c
#include "impcom.inc"
c
      double precision cpu,time
      integer time_array(8)
#if defined (PCWINNT_SYS)
      call cpu_time (cpu)
#else
      real tarray(2),etime
      data tarray /0.0,0.0/
      cpu = dble(etime(tarray))
#endif

      call date_and_time(values=time_array)

      time = time_array(5)*3600 + time_array(6)*60
     &     + time_array(7) + 0.001*time_array(8)
     
      return
      end
