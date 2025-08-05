      subroutine AdjNextTimeStep
c --- primary function: adjust next time increment
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
#include "conv_stat.inc"
#include "dyn_info.inc"
c
      logical  IfEqvDP

      integer j,ieq,inc,iott,neqmax,GetWriteUnit

      double precision
     &  fact2,dtmin,dtmax,dt0,dt1,dt2,dt3,dt4,dt5,dt6,tisave,
     &  con,con1,con2,con3,conts,c1,c2


#if defined(DEBUG)
      call SubBeg('AdjNextTimeStep')
#endif

c --- see if at end of step
      tisave = A_TimeInc
      iott = GetWriteUnit()
      
      if (IfEqvDP(A_TimeCur,A_TimeEnd)) goto 999

c --- set minimun time step
c --- user prescribed
      dtmin = A_MinTimeStep
c --- limit change to factor of 3
      conts = 0.333333333d0*A_TimeInc
      if (conts.gt.dtmin) dtmin = conts

c --- set maximum time step
c --- user prescribed
      dtmax = A_MaxTimeStep
c --- limit change to factor
      conts = 3.0d0*A_TimeInc
      if (conts.lt.dtmax) dtmax = conts
      conts = A_TimeEnd - A_TimeCur
      if (dtmax.gt.conts) dtmax = conts

c --- do not undercut minimum
      if (dtmax.lt.dtmin) dtmin = dtmax

c --- trial time 0: misc. limitations
      dt0 = conts
      if (A_SubStepNum.eq.1) then
c --- do not increase until second substep
          dt0 = A_TimeInc
      elseif (A_TimeInc.lt.A_TimeIno) then
c --- do not increase if trend is decreasing
          dt0 = A_TimeInc
      elseif (A_BisecNum.ge.1) then
c --- do not increase if bisection was active
          dt0 = A_TimeInc
      endif

c --- trial time 1: saturation of equilibrium iterations
      neqmax = A_IterMax
c --- find the max number and max change
      ieq = 1
      inc = 0

      if (ieq.lt.C_MaxIter(1)) ieq = C_MaxIter(1)
      j = C_MaxIter(1) - C_MaxIter(2)
      if (inc.lt.j) inc = j

c --- set min/max time factor
      if (A_Trans.eq.1) then
          fact2 = 2.0d0
      else
          fact2 = 1.50d0
      endif

c --- does not increase by more than 1 iter, open 50%.
      dt1 = A_TimeInc
      if (ieq.le.C_IterCtrlG .and. inc.le.1) dt1 = fact2*A_TimeInc

c --- decrease by 25%.
      if (ieq.gt.C_IterCtrlL) dt1 = 0.75d0*A_TimeInc

c --- trial time 2: time for gap to change status
      dt2 = A_TimeEnd - A_TimeBeg
      if (C_tGap.gt.0.0d0) then
          con = A_TimeEnd - C_tGap
          if (con.lt.dtmin .and. C_StepFactGap.ge.0.9d0) then
              dt2 = A_TimeEnd - A_TimeCur
          else
              dt2 = C_StepFactGap*(C_tGap - A_TimeCur)
          endif
          con = dt2/A_TimeInc
          if (con.le.1.1d0 .and. con.ge.0.9d0) dt2 = A_TimeInc
          write (iott,3000) C_ElGap,C_tGap
 3000     format (' *** ESTIMATED TIME FOR ELEMENT',i8,
     &            ' TO CHANGE STATUS IS ',g11.4)
      endif

c --- trial time 3: negative pivots at the minimum timestep
      dt3 = A_TimeEnd - A_TimeBeg
      if (IfEqvDP(tisave,A_MinTimeStep) .and. C_nPivots.gt.0) then
c --- this is to prevent unnecessary bisecting (we need a better logic here)
          dt3 = A_MinTimeStep
      endif
      
c --- trial time 4: response frequency (2nd order)
      dt4 = A_TimeEnd - A_TimeBeg
      if (A_Trans.eq.1 .and. I_RespFrq.gt.TINY) then
          con = 1.0d0/(I_RespFrq*A_TimeInc)
c --- defaults: steps/cycle = 13 (+10 or -4)
          con1 = 9.0d0
          con2 = 13.0d0
          con3 = 23.0d0
          if (con.lt.con1 .or. con.gt.con3) then
              dt4 = A_TimeInc*con/con2
          else
              dt4 = A_TimeInc
          endif
      endif

c --- time trial 5: creep ratio
      dt5 = A_TimeEnd - A_TimeBeg
      if (A_KeyCreep.ne.0 .and. A_CrMax.gt.TINY
     &    .and. A_CrLimit.gt.0.0d0) then
c --- implicit creep
          con = A_CrLimit/A_CrMax
          c1  = 1.0d0 - A_CrTol
          c2  = 1.0d0 + A_CrTol
c --- only change if changes by more than 10 percent
          if (con.lt.c1 .or. con.gt.c2) then
              dt5 = A_TimeInc*con
          else
              dt5 = A_TimeInc
          endif
      endif

c --- time trial 6: plasticity adjustment
      dt6 = A_TimeEnd - A_TimeBeg
      if (A_KeyPlast.gt.0 .and. A_StraMax.gt.TINY) then
          con = A_StraLimit/A_StraMax
          c1  = 1.0d0 - A_StraTol
          c2  = 1.0d0 + A_StraTol
c --- only change if changes by more than 10 percent
c          write (6,*) 'con, c1, c2', con, c1, c2
          if (con.lt.c1 .or. con.gt.c2) then
              dt6 = A_TimeInc*con
          else
              dt6 = A_TimeInc
          endif
      endif


c --- set minimum time increment
      A_TimeInc = dt0

      if (dt1.lt.A_TimeInc) A_TimeInc = dt1
      if (dt2.lt.A_TimeInc) A_TimeInc = dt2
      if (dt3.lt.A_TimeInc) A_TimeInc = dt3
      if (dt4.lt.A_TimeInc) A_TimeInc = dt4
      if (dt5.lt.A_TimeInc) A_TimeInc = dt5
      if (dt6.lt.A_TimeInc) A_TimeInc = dt6

c --- check lower/upper bounds
      if (A_TimeInc.lt.dtmin) A_TimeInc = dtmin
      if (A_TimeInc.gt.dtmax) A_TimeInc = dtmax

c --- check end of load step and reset (split difference, conts contains time remaining, conts = A_TimeEnd - A_TimeCur)
      con = conts - A_TimeInc
      if (con.lt.dtmin) then
c --- the step is lower than the minimum, so it is OK
          A_TimeInc = conts
      else
c --- OBS: never increase the time-step after reaching the minimum timestep              
          if (abs(A_SpecVer).eq.4) then
              if (A_TimeInc.gt.tisave .and. tisave.eq.dtmin) then
                  A_TimeInc = tisave
              endif
          endif
      endif
c
      if (A_TimeInc.gt.conts) A_TimeInc = conts

c --- print summary
      con = A_TimeInc/tisave
      if (con.ge.1.00001d0) then
          write (iott,3001) A_TimeInc,con
 3001     format (' --> AUTO STEP INCREASED:  NEXT TIME INC =',g12.5,
     &            '  INCREASE FACTOR = ',f7.4)
      elseif (con.le.0.99999d0) then
          write (iott,3002) A_TimeInc,con
 3002     format (' --> AUTO STEP DECREASED:  NEXT TIME INC =',g12.5,
     &            '  DECREASE FACTOR =',f7.4)
      else
          write (iott,3003) A_TimeInc
 3003     format (' --> AUTO STEP UNCHANGED:  NEXT TIME INC =',g12.5)
      endif

 999  continue

#if defined(DEBUG)
      call SubEnd('AdjNextTimeStep')
#endif

      return
      end
