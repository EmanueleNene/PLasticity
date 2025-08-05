      subroutine ContDamage (nvmax,cinem,qa,qb,kerr,presb,taurb,
     &                       tausb,presc0,taurc0,tausc0,rigp0,rigt0,
     &                       idbmode,dfunc,ONED)
c --- primary function: calculate contact pressure and evaluate damage variable
c
c  input arguments:
c       nvmax   (int,sc,in)         - number of state variables
c       cinem   (dp,ar(39),in)      - kinematic parameters and material data
c                                     used for contact and debonding
c       qa      (dp,ar(*),in)       - stresses and state variables at initial
c                                     configuration
c       kerr    (int,sc,in)         - error code
c       presb   (dp,sc,inout)       - normal contact pressure
c       taurb   (dp,sc,inout)       - tangent contact stress in r-dir
c       tausb   (dp,sc,inout)       - tangent contact stress in s-dir
c       rigp0   (dp,sc,in)          - initial normal penalty stiffness
c       rigt0   (dp,sc,in)          - initial tangent penalty stiffness
c       idbmode (int,sc,inout)      - failure mode for debonding
c
c  output arguments:
c       qb      (dp,ar(*),out)      - stresses and state variables at updated
c                                     configuration
c       presc0  (dp,sc,out)         - normal contact pressure w/o damage
c       taurc0  (dp,sc,out)         - tangent contact stress in r-dir w/o damage
c       tausc0  (dp,sc,out)         - tangent contact stress in s-dir w/o damage
c       dfunc   (dp,sc,out)         - damage function value at current substep
c
#include "impcom.inc"
#include "solu_info.inc"
c
      integer nvmax,kerr,idbmode
c
      double precision cinem(41),qa(nvmax),qb(nvmax),presb,taurb,
     &  tausb,presc0,taurc0,tausc0,rigp0,rigt0,dfunc

      integer icflg
c
      double precision ZERO,ONE,tauab,presa,taura,tausa,
     &  presam,tauaam,presc,taurc,tausc,tauac,tauacm,
     &  tnmax,ttmax,gnmax,gtmax,vrpar,
     &  dnew,dold,dvnew,dvold,dfunca,dfuncold,
     &  sfratio,fact1,fact2,vfac1,vfac2,
     &  dtest1,dtest2,ASMALL,ONED
c
      parameter (ZERO=0.0d0,ONE=1.0d0,ASMALL=1.d-8)
      
#if defined(DEBUG)
      call SubBeg('ContDamage')
#endif

c --- initialize
c
      dnew = ZERO
      dold = ZERO
      dvnew = ZERO
      dvold = ZERO

      sfratio = ZERO

      presc = ZERO
      taurc = ZERO
      tausc = ZERO
      presa = ZERO
      taura = ZERO
      tausa = ZERO
      presam = ZERO
      tauaam = ZERO
      tauacm = ZERO
      tauac = ZERO
      tauab = sqrt(taurb*taurb + tausb*tausb)

      dfunc  = ZERO
      dfunca = ZERO
      dfuncold = ZERO

      fact1 = ZERO
      fact2 = ZERO
      vfac1 = ZERO
      vfac2 = ZERO

      dtest1 = ZERO
      dtest2 = ZERO

c --- fracture model data for debonding
c
      tnmax = cinem(36)
      gnmax = cinem(37)
      ttmax = cinem(38)
      gtmax = cinem(39)
      vrpar = cinem(40)
      icflg = nint(cinem(41))

c --- debonding model is based on damage with softening.
c --- debonding process is modelled as initially elastic loading
c --- up to maximum stress (bond strength) followed by linear
c --- softening characterized by a damage variable. The damage
c --- is ZERO at the start of softening and becomes ONE at failure
c --- when the stress goes to zero. Damage is assumed to be
c --- irreversible so if unloading occurs during softening, the
c --- bond stiffness is reduced and when loading goes to zero
c --- both stress and displacement go to zero. Reloading occurs
c --- along the same path as loading (same reduced bond stiffness)
c --- no damage in case of compressive normal stress is assumed 
c --- if icflg.eq.0
c      write (6,*) 'presb ini',presb
      
      if (presb.gt.ZERO) then
          if ((idbmode.eq.1) .or.
     &        (idbmode.eq.3 .and. icflg.eq.0)) then
               dold = qa(nvmax-9)
               dvold = qa(nvmax-11)
               dnew = dold
               dvnew = dvold
               dfunc = ZERO
               presc0 = presb
               taurc0 = taurb
               tausc0 = tausb
               goto 200
          endif
      endif

c --- evaluate damage/failure function based on mode
c
      if (idbmode.eq.1) then
          dfunc = abs(presb)/tnmax - ONE
          sfratio = ONE - tnmax/(rigp0*gnmax)
      elseif (idbmode.eq.2) then
          dfunc = tauab/ttmax - ONE
          sfratio = ONE - ttmax/(rigt0*gtmax)
      elseif (idbmode.eq.3) then
          if (presb.gt.ZERO) then
              fact1 = ZERO
          else
              fact1 = (presb*presb)/(tnmax*tnmax)
          endif
          fact2 = (tauab*tauab)/(ttmax*ttmax)
          dfunc = sqrt(fact1+fact2)-ONE
c --- keep ratio same in both directions by
c --- scaling penalty stiffness for tangential direction
          sfratio = ONE - tnmax/(rigp0*gnmax)
      endif

c --- initialize damage variable for current substep
c     based on failure function value
c
      if (dfunc.lt.ZERO) then
          dfunc = ZERO
          dnew = ZERO
      else
          dnew = dfunc/((dfunc+ONE)*sfratio)
          if (dnew.lt.ZERO) then ! no new damage
              dnew = ZERO
          elseif (dnew.gt.ONED) then ! failure
              dnew = ONE
          endif
      endif

c --- retrieve damage variable from previous substep
c
      dold = qa(nvmax-9)
      dvold = qa(nvmax-11)
      if (dold.gt.ONED) then     ! previously failed
          dnew = ONE
      endif

c --- evaluate visc. regularized damage variable for current substep
c
      vfac1 = A_TimeInc/(A_TimeInc+vrpar)
      vfac2 = vrpar/(A_TimeInc+vrpar)
      dvnew = vfac1*dnew+vfac2*dvold
      if (dvnew.gt.ONED) then
          dvnew = ONE
      endif

c --- update damage variable and modify contact stresses to reflect
c     damage/softening
c
      if ((dvold.eq.ZERO .and. dvnew.gt.ZERO) .or.
     &     dvold.gt.ZERO) then
         presc = presb
         taurc = taurb
         tausc = tausb
         tauac = sqrt(taurc*taurc + tausc*tausc)
         if (dvold.eq.ZERO) then
c --- no previous damage
c
            tauacm = tauab
            if (presc.lt.ZERO) then
                presc = presc*(ONE-dvnew)
            endif
            taurc = taurc*(ONE-dvnew)
            tausc = tausc*(ONE-dvnew)
         elseif (dvold.lt.ONED) then
c --- previously damaged
c
            if (dold.lt.ONED) then
c --- evaluate damage function for last substep with damage
c
               dfuncold = dold/(ONE/sfratio-dold)

c --- evaluate damage function for current substep
c
               taura = qa(2)
               tausa = qa(3)
               tauaam = sqrt(taura*taura + tausa*tausa)/(ONE-dvold)
               tauacm = tauac
               if (idbmode.eq.1) then
                   dfunc = abs(presc)/tnmax - ONE
               elseif (idbmode.eq.2) then
                   dfunc = tauacm/ttmax - ONE
               elseif (idbmode.eq.3) then
                   if (presb.gt.ZERO) then
                       fact1 = ZERO
                   else
                       fact1 = (presc*presc)/(tnmax*tnmax)
                   endif
                   fact2 = (tauacm*tauacm)/(ttmax*ttmax)
                   dfunc = sqrt(fact1+fact2)-ONE
               endif
               if (dfunc.lt.ZERO) then
                   dfunc = ZERO
               endif

c --- evaluate damage function for previous substep
c
               presa = qa(1)
               if (presa.lt.ZERO) then
                   presam = presa/(ONE-dvold)
                   if (presc.lt.ZERO) then
c --- idbmode_old = idbmode_new
                      if (idbmode.eq.1) then
                          dfunca = abs(presam)/tnmax - ONE
                      elseif (idbmode.eq.2) then
                          dfunca = tauaam/ttmax - ONE
                      elseif (idbmode.eq.3) then
                          fact1 = (presam*presam)/
     &                            (tnmax*tnmax)
                          fact2 = (tauaam*tauaam)/
     &                            (ttmax*ttmax)
                          dfunca = sqrt(fact1+fact2)-ONE
                      endif
                  else
c --- idbmode_old = 1, 2 or 3, idbmode_new = 2 or 3
                      if (tnmax.lt.TINY) then
c --- idbmode_old = 2
                          dfunca = tauaam/ttmax - ONE
                      else
                          if (ttmax.lt.TINY) then
c --- idbmode_old = 1
                              dfunca = abs(presam)/tnmax - ONE
                          else
c --- idbmode_old = 3
                              fact1 = (presam*presam)/(tnmax*tnmax)
                              fact2 = (tauaam*tauaam)/(ttmax*ttmax)
                              dfunca = sqrt(fact1+fact2)-ONE
                          endif
                      endif
                  endif
                  if (dfunca.lt.ZERO) then
                      dfunca = ZERO
                  endif
               elseif (presa.gt.ZERO .and.
     &                (idbmode.eq.2 .or.
     &                (idbmode.eq.3 .and. icflg.eq.1))) then
                       dfunca = tauaam/ttmax - ONE
                       if (dfunca.lt.ZERO) then
                           dfunca = ZERO
                       endif
               else
                  dfunca = ZERO
               endif

c --- compare damage functions to determine new damage
c
               dtest1 = abs((dfunca/dfuncold) - ONE)
               if (dtest1.le.ASMALL) then ! damage in pr. subst.
                   if (dfunc.gt.dfunca) then ! new damage/fail.
                       dnew = dfunc/((dfunc+ONE)*sfratio)
                       dtest2 = abs((dnew/dold) - ONE)
                       if (dnew.lt.dold .and. dtest2.gt.ASMALL) then
                           call erhandler ('ContDamage',5000,kerr,
     &                       'Strange - no damage under increased
     &                        stresses and previous damage',0,' ')
                       elseif (dnew.gt.ONED) then ! failure
                           dnew = ONE
                       endif
                  else          !no new damage
                       dnew = dold
                  endif
               elseif (dfunca.lt.dfuncold) then ! no pr. damage
                   if (dfunc.gt.dfunca) then ! loading
                       if (dfunc.lt.dfuncold) then ! no new damage
                           dnew = dold
                       else       !new damage/failure
                           dnew = dfunc/((dfunc+ONE)*sfratio)
                           dtest2 = abs((dnew/dold) - ONE)
                           if (dnew.lt.dold .and.
     &                         dtest2.gt.ASMALL) then
                               call erhandler ('ContDamage',5000,kerr,
     &                         'Strange - no damage under increased
     &                         stresses',0,' ')
                        elseif (dnew.gt.ONED) then ! failure
                           dnew = ONE
                        endif
                     endif
                   else          ! unloading
                      dnew = dold
                   endif
               else
                   dnew = ONE
c --- understand the source of this error now, do the bisection
c                   call PutErr(1)
c                   call erhandler ('ContDamage',5000,3,
c     &                             'Problem with contact stabilization',
c     &                             0,' ')
               endif
            endif

c --- evaluate visc. regularized damage variable
c --- based on update damage variable for current substep
            dvnew = dnew*vfac1+dvold*vfac2
            if (dvnew.gt.ONED) then
                dvnew = ONE
            endif

c --- modify contact stresses to reflect new damage/softening
c
            if (presc.lt.ZERO) then
                presc = presc*(ONE-dvnew)
            endif
            if (tauac.gt.TINY) then
                taurc = taurc*(ONE-dvnew)*(tauacm/tauac)
                tausc = tausc*(ONE-dvnew)*(tauacm/tauac)
            endif
         endif

c --- save contact stresses to SVRs
c
         presc0 = presb
         if (abs(tauab).gt.ASMALL) then
             taurc0 = tauacm*taurb/tauab
             tausc0 = tauacm*tausb/tauab
         else
             taurc0 = ZERO
             tausc0 = ZERO
         endif
         if (dvnew.eq.ONE) then !failure
             if (presc.le.ZERO) then
                 presb = ZERO
             endif
             taurb = ZERO
             tausb = ZERO
         else                   !no/some damage
             presb = presc
             taurb = taurc
             tausb = tausc
         endif
         qb(1) = presb
         qb(2) = taurb
         qb(3) = tausb
      endif
 200  continue

c --- save damage variable for current substep and
c --- current step time to SVRs
c
      if (dvnew.gt.ZERO .and. dvnew.le.ONED) then
          qb(nvmax-8) = A_TimeCur
      endif

      qb(nvmax-9) = dnew
      qb(nvmax-11) = dvnew
      
c      write (6,*) 'presb,presc',presb,presc
     
#if defined(DEBUG)
      call SubEnd('ContDamage')
#endif

      return
      end
