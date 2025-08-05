      subroutine AutoTimeStep (dt,dtmn,dtmx)
c --- primary function: determine the time stepping sizes
c
c  output arguments:
c     dt       (dp,sc,out)         - initial time step size
c     dtmn     (dp,sc,out)         - minimum time step size
c     dtmx     (dp,sc,out)         - maximum time step size
c
#include "impcom.inc"
#include "solu_info.inc"
#include "conv_stat.inc"
c
      logical  IfEqvDP

      integer j
      double precision dt,dtmn,dtmx, timels,con

#if defined(DEBUG)
      call SubBeg('AutoTimeStep')
#endif

c --- define time increments
      timels = A_TimeEnd - A_TimeBeg
      if (IfEqvDP(A_TimeBeg,A_TimeEnd) .or.
     &    A_TimeEnd.lt.A_TimeBeg) then
          timels = 1.0d0
      endif

      dt = timels

c --- set initial time step
      if (C_NSubIni.gt.0) then
c --- check for nsubst input
          con = timels/C_NSubIni
          if (con .lt. dt) dt = con
      endif
c --- reset values for even time stepping
c --- fix for very small time step vs. total time
      con = timels/dt
      if (con .le. 1.0d8) then
          j = con + 0.5d0
          dt = timels/j
      else
          con = mod(timels,dt)
          dt = dt*timels/(timels - con)
      endif

c --- set max/min time step limits if auto time step is used
      if (A_KeyAutoStep.eq.1) then
          dtmn = dt
          dtmx = timels
c --- max time step set
          if (C_NSubMin .gt. 0) then
c --- check for nsubst input
              con = timels/C_NSubMin
              if (con .lt. dtmx) dtmx = con
          endif
c --- min time step set
          if (C_NSubMax .gt. 0) then
c --- check for nsubst input
              con = timels/C_NSubMax
              if (con .lt. dtmn) dtmn = con
          endif
      endif

#if defined(DEBUG)
      call SubEnd('AutoTimeStep')
#endif

      return
      end
