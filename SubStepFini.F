      subroutine SubStepFini(t0,tw0,tc,twc,kreturn)
c --- primary function: finish up a substep
c
#include "impcom.inc"
#include "solu_info.inc"
#include "dyn_info.inc"
#include "RST_file.inc"
c
      integer GetWriteUnit
      integer iott,kreturn
      double precision  timsav, t0, t1, tw0, tw1, tc, twc, con1, con2
      character*8 wall

#if defined(DEBUG)
      call SubBeg('SubStepFini')
#endif

      iott = GetWriteUnit()

c --- write out substep trailer
      call KCL_time(t1,tw1)
c --- increminet current time
      twc = twc + (tw1 - tw0)
c --- correct for the day change
      if ((tw1 - tw0).lt.0) then
          twc = twc + 3600*24
      endif
c --- update previous time
      tw0 = tw1
      write (wall,2000) twc
 2000 format (f8.1)

c --- same as above for CPU time      
      tc = tc + (t1 - t0)
      t0 = t1
      
      write (iott,2010) A_SubStepNum,tc,adjustl(wall),A_CurLoadFact
 2010 format (2x,'-> SUBSTEP',i6,' CONVERGED.',
     &        ' ELPS ',f9.1,' / ',A8,' LOAD FACTOR =',g12.5)

c --- write out transient info
      if (A_Trans.gt.0 .and. I_RespFrq.gt.TINY) then
          con1 = 1.0d0/I_RespFrq
          con2 = con1/A_TimeInc
          write (iott,2020) I_RespFrq,con1,con2
 2020     format (' --> SYSTEM FREQUENCY =',g11.4,
     &            '  PERIOD = ',g11.4,'  PTS/CYC =',g11.4)
      endif

      if (A_MaxEqStra.gt.TINY) then
          write (iott,2021) A_MaxEqStra
 2021     format (' --> Maximum equivalent plastic strain =',g11.4)
      endif

      if (A_CrMax.gt.TINY) then
          write (iott,2022) A_CrMax
 2022     format (' --> Maximum equivalent creep strain =',g11.4)
      endif
      
c --- report current or incremental reaction sum
      if (A_TotReactCur.gt.0.0d0 .and. A_SubStepNum.gt.1) then
          if (R_OutCode.gt.0 .and. A_SubStepNum.gt.2) then
              con1 = A_TotReactCur - A_TotReactMax
              write (iott,2030) con1, con1/A_TotReactMax*100.0d0
 2030         format (' --> REACTION INC =',g11.3,' (',f6.1,'%)')
c --- force output during next step and exit solution
              if (A_Fail.gt.0) then
c --- allowable difference
                  con2 = A_Fail/100.0d0
c --- actual difference
                  con1 = con1 / A_TotReactMax
                  if (con1.lt.0.0d0 .and. abs(con1).gt.con2 .and. 
     &                A_CurLoadFact.gt.0.1d0) then
c --- failure point reached - quit
                      R_Stop = -2
                      kreturn = 999
                  endif              
              endif
          else
              write (iott,2040) A_TotReactCur
 2040         format (' --> CURRENT REACTION =',g11.3)
          endif 
          A_TotReactMax = max(A_TotReactMax,A_TotReactCur)
      endif

      write (iott,'(/)')
      
      if (kreturn.eq.999) goto 999

c --- nonlinear stabilization information
      A_StabFirstStep = 0

c --- auto load adjustment to time step
      timsav = A_TimeInc

      if (A_KeyAutoStep.eq.1) call AdjNextTimeStep

      A_TimeIno = timsav
      
 999  continue

#if defined(DEBUG)
      call SubEnd('SubStepFini')
#endif

      return
      end      