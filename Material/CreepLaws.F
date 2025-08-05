      subroutine CreepLaws (cropt, propcr , time, dtime,
     &                      moist, MoistOffst, creq_t,
     &                      dpcr, dcrds, creq, sigy)
c --- primary function: calculate incremental creep strain, derivitives of creep
c                       strain based on defined creep law
c
c  input arguments:
c      time       (dp,sc,in)             - current time
c      dtime      (dp,sc,in)             - current time increment
c      propcr     (dp,ar(*),in)          - mat properties array
c      cropt      (in,sc,in)             - creep model index
c      sigy       (dp,sc,in)             - equivalent effective stress at 't' / 't+dt'
c      moist      (dp,sc,in)             - current moisture increment
c      MoistOffst (dp,sc,in)             - moisture offset
c      creq_t     (dp,sc,in)             - equivalent effective creep strain at 't'
c      creq       (dp,sc,in)             - equivalent effective creep strain at 't+dt'
c
c  output arguments:
c      dpcr       (dp,sc,o)              - incremental creep strain
c      dcrds      (dp,ar(*),o)           - output array:
c                                          dcrds(1) - derivative of incremental creep 
c                                                     strain to effective stress
c                                          dcrds(3) - derivative of incremental creep 
c                                                     strain to creep strain
c
#include "impcom.inc"
#include "constants.inc"
c

      double precision ONE, TWO,
     &                 ZERO, THREE,
     &                 ln10
      parameter       (ZERO  = 0.d0,
     &                 ONE   = 1.d0,
     &                 TWO   = 2.d0,
     &                 THREE = 3.d0,
     &                 ln10 = 2.302585093d0)
c

c

      integer          cropt
      double precision dtime,  moist, MoistOffst,
     &                 creq, dpcr, sigy, creq_t
      double precision propcr(*), dcrds(*)
c

c

      double precision time,
     &                 c1, c2, c3, c4, c5, c6,
     &                 c7, c8, c9, c10, c11, c12,
     &                 con1, con2, con3, con4, con5, con6,
     &                 con7, con8, r, drds, T, sigy1,
     &                 logsy, sqTiny, expt
c

      sqTiny   = sqrt (TINY)
c
      dpcr = ZERO
      call v_zero(dcrds(1),4)
      if (sigy. LE. sqTiny .and. creq .le. sqTiny) then
         sigy   = ZERO
         creq   = ZERO
         goto 999
      endif

      if (dtime .le. 0.0d0) then
c --- No creep
          goto 999
      elseif (propcr(1) .eq. 0.0d0 .and. cropt .lt. 13) then
c --- No creep is considered
          goto 999
      elseif (cropt .le. 0) then
c --- No creep is considered
          goto 999
      endif

      T = moist + MoistOffst
c
      goto ( 10,  20,  30,  40,  50,  60,  70,  80,  90, 100,
     &      110, 120, 130, 140, 150), cropt
   10 continue
c --- Primary creep function  C6 = 0
c     dotcreq := c1 * q ^ n * creq ^ m * exp (-b/T)
      c1 = propcr(1)
      c2 = propcr(2)
      c3 = propcr(3)
      c4 = propcr(4)
      expt = ONE
      if (abs(c4).gt.TINY.and.abs(T).gt.TINY) then
          expt = exp(-c4/T)
      endif
      if (creq .le. ZERO) then
          if (c3 .ne. ONE) then
              dpcr = exp((log(c1) + c2 * log(sigy) + log(dtime))/
     &               (ONE - c3))
          else
             dpcr = exp(exp((log(c1) + c2 * log(sigy) +
     &              log(dtime))))
          endif
          dpcr = dpcr * expt
          dcrds(1) = c2 * dpcr / sigy
          dcrds(3) = ZERO
      else
         if (c1.gt.ZERO) dpcr = (exp( log(c1) +  c2 * log(sigy) +
     &                           c3 * log(creq) )) * expt * dtime
         dcrds(1) = c2 * dpcr / sigy
         dcrds(3) = c3 * dpcr / creq
      endif
      goto 999
   20 continue
c --- Primary creep function  C6 = 1
c     dotcreq := c1 * q ^ n * t ^ m * exp (-b/T)
      c1 = propcr(1)
      c2 = propcr(2)
      c3 = propcr(3)
      c4 = propcr(4)
      expt = ONE
      if (abs(c4).gt.TINY.and.abs(T).gt.TINY) then
          expt = exp(-c4/T)
      endif
      if (c1.gt.ZERO) dpcr = (exp ( log(c1) + c2 * log(sigy) +
     &                        c3 * log(time) )) * expt * dtime
      dcrds(1) = c2 * dpcr / sigy
      goto 999
   30 continue
c --- Primary creep function  C6 = 2
c     dotcreq := c1 * q ^ n * r * t ^ m * exp (-r*t)
c     r       := c5 q ^ c3 * exp(-b/T)
      c1 = propcr(1)
      c2 = propcr(2)
      c3 = propcr(3)
      c4 = propcr(4)
      c5 = propcr(5)
      expt = ONE
      if (abs(c4).gt.TINY.and.abs(T).gt.TINY) then
          expt = exp(-c4/T)
      endif
      sigy1   = ONE / sigy
      if (c1 .gt. ZERO.and.c5 .gt. ZERO) then
          r = exp( log(c5) + c3 * log(sigy)) * expt
          drds = c3 * r * sigy1
          dpcr = (exp(log(c1) + c2 * log(sigy) + log(r)
     &         - r * time)) * dtime
          dcrds(1) = (c2 * sigy1 + (ONE /r  - time) * drds ) * dpcr
      endif
      goto 999
   40 continue
c --- 8 parameters creep law
c --- dotepcr = c1 * sigy^c2 * (t^c3+c4 * t^c5 + c6 * t^c7) exp(-c8/T)
      c1 = propcr(1)
      c2 = propcr(2)
      c3 = propcr(3)
      c4 = propcr(4)
      c5 = propcr(5)
      c6 = propcr(6)
      c7 = propcr(7)
      c8 = propcr(8)
      expt = ONE
      if (abs(c8).gt.TINY.and.abs(T).gt.TINY) then
          expt = exp(-c8/T)
      endif
      if (c1.gt.ZERO) dpcr =  exp( log(c1) + c2 * log(sigy) +
     &           log(time ** c3 + c4 * time ** c5 + c6 * time ** c7))
      dpcr = dpcr * expt * dtime
      dcrds(1) = c2 * dpcr / sigy
      goto 999
   50 continue
c --- 7 parameters creep law
      c1 = propcr(1)
      c2 = propcr(2)
      c3 = propcr(3)
      c4 = propcr(4)
      c5 = propcr(5)
      c6 = propcr(6)
      c7 = propcr(7)
      con1 = ZERO
      con2 = ZERO
      con3 = ZERO
      if (c1.gt.ZERO) con1 = exp (log(c1) + c2 * sigy)
      if (c3.gt.ZERO) con2 = exp (log(c3) + c5 * log(sigy/c4))
      if (c6.gt.ZERO) con3 = exp (log(c6) + c7 * sigy)
      con4 = exp ( - con2 * time )
      dpcr = (con1 * ( ONE - con4 ) + con3 * time) * dtime
      con5 = c2 * con1
      con6 = c5 * con2 / sigy
      con7 = c7 * con3
      dcrds(1) = (con5 * (ONE - con4) + con1 * time * con4 * con6 +
     &            time * con7) * dtime
      goto 999
   60 continue

c --- Time hardening creep law for constant stress only
c     This is identical to cropt = 2, except assuming
c     constant stress so that creep strain is explicit
c     integrated.
c     dotepcr := c1 * q ^ n * t ^ m
c     epcr    := c1 * q ^ n * t ^ (m+1) / (1 + m)
c
      c1 = propcr(1)
      c2 = propcr(2)
      c3 = propcr(3)
      c4 = propcr(4)
      expt = ONE
      if (abs(c4).gt.TINY.and.abs(T).gt.TINY) then
          expt = exp(-c4/T)
      endif
      if (c1.gt.ZERO)
     &    dpcr   = exp(log(c1)+c2*log(sigy)) * (time ** (ONE+c3)
     &           - (time - dtime) ** (ONE+c3)) /(ONE+c3)
      dpcr = dpcr * expt
      dcrds(1) = c2 * dpcr / sigy
      goto 999
   70 continue
c     dotepcr := {c1 q^n [(m+1) epcr ]^m} 1/(m+1)
      c1 = propcr(1)
      c2 = propcr(2)
      c3 = propcr(3)
      c4 = propcr(4)
      expt = ONE
      if (abs(c4).gt.TINY.and.abs(T).gt.TINY) then
          expt = exp(-c4/T)
      endif
      if (creq .le. ZERO) then
          creq = exp(log(c1)+c2 * log(sigy) +
     &          (ONE+c3)*log(time))/(ONE+c3)
      endif
      if (c1.gt.ZERO) con1 = (log(c1) + c2 * log(sigy) +
     &                        c3 * log((ONE+c3) * creq))/(ONE+c3)
      dpcr = exp(con1) * expt * dtime
      con1 = dpcr / (c3 + ONE)
      dcrds(1) = c2 * con1 / sigy
      dcrds(3) = c3 * con1 / creq
      goto 999
   80 continue
c --- Hyperbolic creep law
c     dotepcr := c1 * [sinh(c3 * sigy)]^c2 exp(-c4/T)
      c1 = propcr(1)
      c3 = propcr(2)
      c2 = propcr(3)
      c4 = propcr(4)
      expt = ONE
      if (abs(c4).gt.TINY.and.abs(T).gt.TINY) then
         expt = exp(-c4/T)
      endif
      con1 = sinh( c3 * sigy )
      con2 = tanh( c3 * sigy )
      if (c1.gt.ZERO)
     &    dpcr = (exp ( log(c1) + c2 * log(con1))) * expt * dtime
      dcrds(1) = c2 * c3 * dpcr / con2
      goto 999
   90 continue
c --- Secondary creep
c      dotepcr := c1 * exp(sigy/c2)*exp(-c3/T)
      c1 = propcr(1)
      c2 = propcr(2)
      c3 = propcr(3)
      expt = ONE
      if (abs(c3).gt.TINY.and.abs(T).gt.TINY) then
          expt = exp(-c3/T)
      endif
      if (c1.gt.ZERO) dpcr = exp(log(c1)+sigy/c2) * expt * dtime
      dcrds(1) = dpcr/c2
      goto 999
  100 continue
c --- Secondary creep
c      dotepcr := c1 * sigy^c2*exp(-c3/T)
      c1 = propcr(1)
      c2 = propcr(2)
      c3 = propcr(3)
      expt = ONE
      if (abs(c3).gt.TINY.and.abs(T).gt.TINY) then
          expt = exp(-c3/T)
      endif
      if (c1.gt.ZERO) dpcr = exp(log(c1)+c2*log(sigy)) * expt * dtime
      dcrds(1) = c2*dpcr/sigy
      goto 999
  110 continue
c --- Bailey and Secondary creep
c      epcr := (c1 * sigy^c2 * t ^c3 + c4 * sigy^c5 * t) * exp(-c6/T)
      c1 = propcr(1)
      c2 = propcr(2)
      c3 = propcr(3)
      c4 = propcr(4)
      c5 = propcr(5)
      c6 = propcr(6)
      c7 = propcr(7)
      expt = ONE
      if (abs(c4).gt.TINY.and.abs(T).gt.TINY) then
          expt = exp(-c4/T)
      endif
      con1 = ZERO
      if (c1.gt.ZERO)
     &    con1 = exp(log(c1)+c2*log(sigy)) * (time ** (ONE+c3)
     &         - (time - dtime) ** (ONE+c3)) /(ONE+c3) * expt
      expt = ONE
      if (abs(c7).gt.TINY.and.abs(T).gt.TINY) then
          expt = exp(-c7/T)
      endif
      con2 = ZERO
      if (c5.gt.ZERO)
     &    con2 = exp(log(c5)+c6*log(sigy)) * dtime * expt
      dpcr = con1 + con2
      dcrds(1) = (c2*con1+c6*con2)/sigy
      goto 999
  120 continue
c --- general material rational polynomial
c      dotepcr := c1(deps_c/dt)*Dt
      c1 = propcr(1)
      c2 = propcr(2)
      if (abs(c2) .gt. TINY) then
          c3 = propcr(3)
          c4 = propcr(4)
          c5 = propcr(5)
          c6 = propcr(6)
          c7 = propcr(7)
          c8 = propcr(8)
          c9 = propcr(9)
          c10  = propcr(10)
          c11  = propcr(11)
          c12  = propcr(12)
          sigy1   = ONE / sigy
          logsy   = log(sigy)
          con1 = exp(log(c2)+c3*sigy*ln10+c4*logsy)
          con2 = exp(log(c7 )+c8 *log(con1)+c9 *logsy)
          con3 = exp(log(c10)+c11*log(con1)+c12*logsy)
          con4 = ONE/(ONE+con3*time)
          con5 = con4*con4
          dpcr = c1*(con2*con3*con5+con1)*dtime
          con6 = (c3*ln10+c4*sigy1)*con1
          con7 = (c8 /con1*con6+c9 *sigy1)*con2
          con8 = (c11/con1*con6+c12*sigy1)*con3
          dcrds(1)= c1*(con7*con3*con5+con2*con4*con5*(ONE-con3*time)*
     &              con8+con6)*dtime
      else
          dpcr = ZERO
          dcrds(1) = ZERO
      endif
      goto 999
  130 continue

c --- Findley creep equation for plastics
      c2 = propcr(1)
      c3 = propcr(2)
      c4 = propcr(3)
      c5 = propcr(4)
      c6 = propcr(5)
      c7 = propcr(6)
      con1 = c2 * sigy + c3 * sigy * sigy + c4 * sigy ** THREE
      con2 = c2 + TWO * c3 * sigy + THREE * c4 * sigy * sigy
      con3 = c5 + c6 * sigy
      expt = ONE
      if (abs(c7).gt.TINY.and.abs(T).gt.TINY) then
          expt = exp(-c7/T)
      endif
      con4 = time ** con3
      dpcr = con1 * con4 - creq_t
      con5 = c6 * log(time)
c     con4 = con3 * con4 / time * dtime
      dcrds(1) = (con2 + con1 * con5) * expt * con4
      goto 999
  140 continue

c --- power law with tertiary creep
c
      c6 = propcr(6)
      expt = ONE
      if (abs(c6).gt.TINY.and.abs(T).gt.TINY) then
          expt = exp(-c6/T)
      endif

c --- power law
      dcrds(1)= ZERO
      c1 = propcr(1)
      c2 = propcr(2)
      con1 = ZERO
      if (c1.gt.ZERO) then
          con1 = exp(log(c1)+c2*log(sigy)) * dtime
          dcrds(1)= c2*con1/sigy
      endif

c --- tertiary creep
      c3 = propcr(3)
      c4 = propcr(4)
      c5 = propcr(5)
      con3 = ZERO
      if(c3.gt.ZERO) then
         con2 = ONE/(ONE-c5)
         con3 = exp(con2*(log(c3)+log(ONE-c5)
     &        + c5*log(time))
     &        + c4*log(sigy))
         con3 = con2 * con3 * dtime
         dcrds(1) = dcrds(1) + c4*con3/sigy
      endif
      dpcr = (con1 + con3) * expt
      dcrds(1) = dcrds(1) * expt

      goto 999
  150 continue

c --- power law with tertiary creep
c
      c6 = propcr(6)
      expt = ONE
      if(abs(c6).gt.TINY.and.abs(T).gt.TINY) then
         expt = exp(-c6/T)
      endif
c
c --- power law
      c1 = propcr(1)
      c2 = propcr(2)
      con1 = ZERO
      if(c1.gt.ZERO) con1 = exp(log(c1)+c2*log(sigy)) * time
      dcrds(1)= c2*con1/sigy

c --- tertiary creep
      c3 = propcr(3)
      c4 = propcr(4)
      c5 = propcr(5)
      con3 = ZERO
      if(c3.gt.ZERO) then
         con2 = ONE/(ONE-c5)
         con3 = exp(con2*(log(c3)+log(ONE-c5)+log(time))
     &           + c4*log(sigy))
      endif
      dpcr = (con1 + con3) * expt - creq_t
      dcrds(1) = dcrds(1) + c4*con3/sigy
      dcrds(1) = dcrds(1) * expt

  999 continue

      return
      end

