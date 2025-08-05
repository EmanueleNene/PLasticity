      subroutine GapTiming (elem,u,ptarget,sfact)
c --- primary function: predict time for gap status changes
c
c  input arguments:
c     elem     (int,sc,in)          - element number
c     u        (dp,ar(3),in)        - values vs time used in the prediction
c     ptarget  (dp,sc,in)           - predicted ptarget value
c     sfact    (dp,sc,in)           - stepping factor (usually 0.5)
c
c  referenced variables:
c     C_tGap        (dp,sc,auto)     - critical gap time prediction
c     C_ElGap       (int,sc,auto)    - element associated with prediction
c     C_StepFactGap (dp,sc,auto)     - auto load step factor
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
#include "conv_stat.inc"
c
      integer elem
c
      double precision
     & u(3),ptarget,sfact,valtim,tmold1,tmold2,valtm1,valtm2,uavg,
     & con1,con2,con3,acon,bcon,ccon
c
#if defined(DEBUG)
      call SubBeg('GapTiming')
#endif

c --- fit gap values versus time
      valtim = -999.0d0
      tmold1 = A_TimeCur - A_TimeInc
      tmold2 = tmold1 - A_TimeIno
      if (A_TimeInc.le.0.0d0 .or. tmold1.le.0.0d0) goto 999
      uavg = 0.33333333d0*(abs(u(1)) + abs(u(2)) + abs(u(3)))
      if (tmold2.gt.0.0d0 .and. A_TimeIno.gt.0.0d0) then

c --- quadratic fit
c
         con1 = u(1)/(A_TimeInc*(A_TimeCur - tmold2))
         con2 = -u(2)/(A_TimeInc*A_TimeIno)
         con3 = -u(3)/((tmold2 - A_TimeCur)*A_TimeIno)
         acon = con1 + con2 + con3
         bcon = -con1*(tmold1+tmold2) - con2*(A_TimeCur+tmold2)
     &          -con3*(A_TimeCur+tmold1)
         ccon = con1*tmold1*tmold2 + con2*A_TimeCur*tmold2
     &        + con3*A_TimeCur*tmold1 - ptarget
         con1 = bcon*bcon - 4.0d0*acon*ccon
c --- check for imaginary roots
         if (con1 .le. 0.0d0) goto 999
         con2 = acon*A_TimeCur*A_TimeCur
         if (uavg .gt. 0.0d0) con2 = con2/uavg
         if (abs (con2) .gt. 1.0d-9) then
c --- quadratic extrapolation
             valtm1 = (-bcon - sqrt(con1))/(2.0d0*acon)
             valtm2 = (-bcon + sqrt(con1))/(2.0d0*acon)
c --- make sure valtm2 > valtm1
             if (valtm2 .lt. valtm1) then
                 con3 = valtm1
                 valtm1 = valtm2
                 valtm2 = con3
             endif
             con1 = A_TimeCur - A_TimeInc/2.0d0
             if (valtm1 .gt. con1) then
                 valtim = valtm1
                 if (valtim .lt. A_TimeCur) valtim = A_TimeCur
             elseif (valtm2 .gt. con1) then
                 valtim = valtm2
                 if (valtim .lt. A_TimeCur) valtim = A_TimeCur
             endif
         else
c --- linear extrapolation
             con2 = bcon*A_TimeCur
             if (uavg .gt. 0.0d0) con2 = con2/uavg
             if (abs (con2) .gt. 1.0d-9) then
                 valtm1 = -ccon/bcon
                 con1 = A_TimeCur - A_TimeInc/2.0d0
                 if (valtm1 .gt. con1) then
                     valtim = valtm1
                     if (valtim .lt. A_TimeCur) valtim = A_TimeCur
                 endif
             endif
         endif
      else
c --- linear fit
         bcon = (u(1) - u(2))/A_TimeInc
         ccon = (u(2)*A_TimeCur - u(1)*tmold1)/A_TimeInc - ptarget
         con2 = bcon*A_TimeCur
         if (uavg .gt. 0.0d0) con2 = con2/uavg
         if (abs (con2) .gt. 1.0d-9) then
c --- linear extrapolation
             valtm1 = -ccon/bcon
             con1 = A_TimeCur - A_TimeInc/2.0d0
             if (valtm1 .gt. con1) then
                 valtim = valtm1
                 if (valtim .lt. A_TimeCur) valtim = A_TimeCur
             endif
         endif
      endif

c --- attempt to set the critical time
      if (valtim .gt. 0.0d0) then
c
          con1 = 1000.0d0*(A_TimeEnd-A_TimeBeg)

          if (valtim .gt. con1) goto 999
          if (C_tGap.le.0.0d0 .or. valtim.lt.C_tGap) then
c --- set the critical time
              call lockSMP (W_LockAutoStep)
              C_tGap = valtim
              C_ElGap = elem
              C_StepFactGap = sfact
              call unlockSMP (W_LockAutoStep)
         endif
      endif

 999  continue

#if defined(DEBUG)
      call SubEnd('GapTiming')
#endif

      return
      end
