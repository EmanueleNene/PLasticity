      subroutine ContPredi (elem,kpred,nvmax,npi,qa,qb,mu,ipair,dmax)
c --- primary function: service routine for contact elements - predicts contact status
c
c  input arguments:
c    elem     (int,sc,in)      - element number
c    nvmax    (int,sc,in)      - number of state variables on each i.p.
c    npi      (int,sc,in)      - number of intergration points
c    qa       (dp,ar(nvmax,npi),in) - stresses and state variables at initial
c                                     configuration for each i.p.
c    qb       (dp,ar(nvmax,npi),in) - stresses and state variables at updated
c                                     configuration for each i.p.
c    mu       (dp,sc,in)       - friction coefficient
c    kpred    (int,sc,in)      - prediction control key (1 or 2)
c    dmax     (dp,sc,in)       - pinball
c    ipair    (int,sc,in)      - contact element numbering
c
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
#include "conv_stat.inc"
c
      integer  elem,ipi,npi,kstat,kstat1,kstat2,kpred,nvmax,ipair
c
      double precision
     &  gap,gap1,gap2,fn,fn1,fn2,velr,velr1,velr2,taur,taur1,taur2,
     &  vels,vels1,vels2,taus,taus1,taus2,dmax,
     &  vect(3),qa(nvmax,npi),qb(nvmax,npi),mu,fact,valtim
c

#if defined(DEBUG)
      call SubBeg('ContPredi')
#endif

      do 100 ipi = 1,npi
         kstat  = nint(qb(5,ipi))
         kstat1 = nint(qa(5,ipi))
         kstat2 = nint(qa(nvmax-6,ipi))
         if (kstat.eq.0) kstat = 4
         if (kstat1.eq.0) kstat1 = 4
         if (kstat2.eq.0) kstat2 = 4
c --- skip out of pinball
         if (kstat.eq.4 .and. kstat1.eq.4) goto 100
         gap  = qb(9,ipi)
         gap1 = qa(9,ipi)
         gap2 = qa(nvmax-5,ipi)
         fn  = qb(1,ipi)
         fn1 = qa(1,ipi)
         fn2 = qb(nvmax-4,ipi)
         velr  = qb(12,ipi)
         velr1 = qa(12,ipi)
         velr2 = qa(nvmax-3,ipi)
         vels  = qb(13,ipi)
         vels1 = qa(13,ipi)
         vels2 = qa(nvmax-2,ipi)
         taur  = qb(2,ipi)
         taur1 = qa(2,ipi)
         taur2 = qa(nvmax-1,ipi)
         taus  = qb(3,ipi)
         taus1 = qa(3,ipi)
         taus2 = qa(nvmax,ipi)
         if (kpred .eq. 2) then
             fact = 0.5d0
         else
             if (kstat .le. 3) then
                 fact = 1.0d0
             else
                 fact = 0.5d0
             endif
         endif
c --- gap to close
         if (kstat.le.2 .and. kstat1.eq.4) then
c --- predict time for contact
             valtim = A_TimeCur + 0.5d0*A_TimeInc
             if (C_tGap.le.0.0d0 .or. valtim.lt.C_tGap) then
                 call lockSMP (W_LockAutoStep)
                 C_tGap = valtim
                 C_ElGap = elem
                 C_StepFactGap = fact
                 call unlockSMP (W_LockAutoStep)
             endif
             goto 100
         endif
         if (kstat.eq.3 .and. kstat1.eq.4.and.ipair.eq.0) then
c --- predict time for contact
             valtim = A_TimeCur + 0.5d0*A_TimeInc
             if (C_tGap.le.0.0d0 .or. valtim.lt.C_tGap) then
                 call lockSMP (W_LockAutoStep)
                 C_tGap = valtim
                 C_ElGap = elem
                 C_StepFactGap = fact
                 call unlockSMP (W_LockAutoStep)
              endif
              goto 100
         endif
         if (kstat.ge.3 .and. kstat1.ge.3) then
c --- predict time for contact
             if (abs(gap).le.1.d-6*dmax) goto 100
             if (gap.lt.gap1 .and. gap1.lt.gap2) then
                 vect(1) = gap
                 vect(2) = gap1
                 vect(3) = gap2
                 call GapTiming (elem,vect(1),0.0d0,fact)
             endif
             goto 100
         endif
c --- gap to open
         if (kstat.eq.4 .and. kstat1.le.2.and.ipair.eq.0) then
             valtim = A_TimeCur + 0.5d0*A_TimeInc
             if (C_tGap.le.0.0d0 .or. valtim.lt.C_tGap) then
                 call lockSMP (W_LockAutoStep)
                 C_tGap = valtim
                 C_ElGap = elem
                 C_StepFactGap = fact
                 call unlockSMP (W_LockAutoStep)
             endif
             goto 100
         endif
         if (kstat.le.2 .and. kstat1.le.2) then
c --- predict time for open
             if (abs(gap).le.1.d-6*dmax) goto 100
             if (-gap.lt.-gap1 .and. -gap1.lt.-gap2) then
                 vect(1) = -gap
                 vect(2) = -gap1
                 vect(3) = -gap2
                 call GapTiming (elem,vect(1),0.0d0,fact)
             endif
             if (fn.lt.fn1 .and. fn1.lt.fn2) then
c --- predict time to open
                 vect(1) = fn
                 vect(2) = fn1
                 vect(3) = fn2
                 call GapTiming (elem,vect(1),0.0d0,fact)
             endif
             if (fn.gt.fn1 .and. fn1.gt.fn2.and.kpred.eq.2) then
c --- predict time to contact
                 vect(1) = fn
                 vect(2) = fn1
                 vect(3) = fn2
                 call GapTiming (elem,vect(1),1.0d0,fact)
             endif
             if (kpred.eq.2 .and. mu .ne. 0.d0) then
c --- sliding occur
                 if (abs(velr)+abs(vels).le.1.d-5*dmax) goto 100
                 if (abs(velr).lt.abs(velr1).and.abs(velr1).lt.
     &                            abs(velr2)) then
c --- predict time to stop sliding
                    vect(1) = abs(velr)
                    vect(2) = abs(velr1)
                    vect(3) = abs(velr2)
                    call GapTiming (elem,vect(1),0.0d0,fact)
                 endif
                 if (abs(vels).lt.abs(vels1).and.abs(vels1).lt.
     &                            abs(vels2)) then
c --- predict time to stop sliding
                    vect(1) = abs(vels)
                    vect(2) = abs(vels1)
                    vect(3) = abs(vels2)
                    call GapTiming (elem,vect(1),0.0d0,fact)
                 endif
                 if (abs(taur).gt.abs(taur1).and.abs(taur1)
     &                        .gt.abs(taur2)) then
c --- predict time to stop sliding
                    vect(1) = abs(taur)
                    vect(2) = abs(taur1)
                    vect(3) = abs(taur2)
                    call GapTiming (elem,vect(1),1.0d0,fact)
                 endif
                 if (abs(taus).gt.abs(taus1).and.abs(taus1)
     &                        .gt.abs(taus2)) then
c --- predict time to stop sliding
                    vect(1) = abs(taus)
                    vect(2) = abs(taus1)
                    vect(3) = abs(taus2)
                    call GapTiming (elem,vect(1),1.0d0,fact)
                endif
             endif
         endif
 100  continue
c

#if defined(DEBUG)
      call SubEnd('ContPredi')
#endif

      return
      end
