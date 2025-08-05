      subroutine CreepPlastDriver (keycut, moist, nTens, elemId, matId,
     &                             prop, propcr, eEl, ePl, eCr,
     &                             eqpl, eqcr, matmisc,
     &                             wt1, wt2, wt3, etaEl,
     &                             stress0, Zee, sigElp,
     &                             stress, dsdePl, nlkword,
     &                             sedEl, sedPl, sedCr, 
     &                             sdv)
c --- primary function: driver for creep and isotropic viscoplasticity for beam elements
c
c --- Reference
c          Based on 'Elastoplasticity and viscoplasticity'
c          Simo, J.C. and Hughes, T.J.R., 1998
c
c  input arguments:
c      moist     (dp,sc,in)                - moisture at the matl point
c      sigElp   (dp,ar(nTens),in)          - stresses from elastic predictor
c      statev0  (dp,ar(*),in)              - internal state variables 't'
c      nlkword (int*8,ar(2,*), out)        - nonlinear material index
c             nlkword(1,1)                 - virtual starting position of plasticity
c             nlkword(2,1)                 - plasticity option:
c             nlkword(1,2)                 - virtual starting position of the data
c             nlkword(2,2)                 - creep option
c
c  output arguments:
c      eEl      (dp,ar(nTens),inout)       -  eEl_t + eIncTot/ eEl_tDt
c      ePl      (dp,ar(nTens),inout)       -  ePl_t/ ePl_tDt
c      stress   (dp,ar(nTesn),inout)       -  stress(elprd) / stress_tDt
c      dsdePl   (dp,ar(nTens,nTens),inout) -  el stiff / el-pl stiff
c      eqpl     (dp,ar(*),inout)           -  eqpl_t / eqpl_tDt
c                                             plWork_t / plWork_tDt
c      wt1      (dp,ar(*),work)            -  work tensor - Sel
c      wt2      (dp,ar(*),work)            -  work tensor - inc pl strain
c      wt3      (dp,ar(*),work)            -  work tensor - inc creep strain
c      Zee      (dp,ar(nTens,*),work)      -  work tensor - 4th order I'
c
c  local variables:
c      sigy_t   (dp,sc,loc)                -  yield stress at 't'
c      gmu      (dp,sc,loc)                -  shear moduli
c      eEldev   (dp,ar(6),loc)             -  deviatoric elastic strain at t+Dt
c      etaEl    (dp,ar(6),loc)             -  deviatoric back-stress    at t+Dt
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"

      double precision HALF,    THIRD,
     &                 ZERO,    THREE, ONEDM12
      parameter       (ZERO       = 0.d0,
     &                 HALF       = 0.5d0,
     &                 THIRD      = 1.d0/3.d0,
     &                 THREE      = 3.d0,
     &                 ONEDM12    = 1.d-12)
c
      double precision EqvStrain

      integer*8  nlkword(W_nMatCol,W_nMatRow)
      integer    nTens, elemId, matId, keycut

      double precision
     &                 dsdePl(nTens,*), propcr(*),
     &                 eEl(*), ePl(*), eCr(*), eqpl(*), eqcr(*),
     &                 sigElp(nTens), etaEl(nTens),  sdv(*),
     &                 wt1(nTens), wt2(nTens), wt3(nTens),
     &                 stress0(nTens), Zee(nTens,nTens),
     &                 stress(*), sedEl, sedPl, sedCr, matmisc(*),
     &                 Flstress(6), moist

      integer          i, cropt, plstate
      integer*8        keyplL, keycrL
      double precision plWork_t, crWork_t, sigy, dperr(3), prop(13),
     &                 creq, dpcr, pleq, dpleq, C(5),
     &                 dtime, time


      keyplL  = nlkword(1,1)
      keycrL  = nlkword(1,2)
      cropt   = nlkword(2,2)
      dtime   = A_TimeInc
      time    = A_TimeCur
      plstate = 0
      dpleq   = ZERO
      dpcr    = ZERO
      pleq    = ZERO
      creq    = ZERO
      sigy    = ZERO
      if (keyplL.gt.0) then
          pleq     = abs(eqpl(1))
          plWork_t = eqpl(2)
      endif
      if (keycrL.gt.0) then
          creq     = eqcr(1)
          crWork_t = eqcr(2)
          sigy     = matmisc(2)
      endif

      C(1) = ZERO
      C(2) = ZERO
      do i = 1, nTens
         C(1) = C(1) + abs( eEl(i) )
         C(2) = C(2) + abs( sigElp(i) )
      enddo

c --- zero out strains
      call v_zero (wt1(1),nTens)
      call v_zero (wt2(1),nTens)
      call v_zero (wt3(1),nTens)
c
      if (C(1)+C(2) .le. ONEDM12) goto 500

      call v_move (stress(1), stress0(1), nTens)

      sigy = sqrt ( stress0(1) * stress0(1)
     &     +  THREE * ( stress0(2) * stress0(2) +
     &                  stress0(3) * stress0(3) ) )
c
      call v_move (sigElp(1), etaEl(1), nTens)

      if (keyplL.gt.0 .and. keycrL.gt.0) then
c --- strain-rate independent plasticity plus creep
          call CrPlMatCal (keycut, elemId, matId,
     &                     prop(1), propcr(1), time, dtime,
     &                     nTens, nlkword(1,1),dsdePl(1,1),etaEl(1),
     &                     wt1(1), wt2(1), wt3(1), Zee(1,1),
     &                     pleq, creq, dpleq, dpcr,
     &                     sigy, moist, A_MoistOffst, plstate)

      elseif (keycrL.gt.0) then
c --- pure creep
          call CrMatCal (keycut,
     &                   prop(1), propcr(1), time, dtime,
     &                   nTens, nlkword(1,1),dsdePl(1,1),etaEl(1),
     &                   wt1(1), wt3(1), Zee(1,1),
     &                   creq, dpcr,
     &                   sigy, moist, A_MoistOffst, plstate)
      elseif (keyplL.gt.0) then
c --- pure plasticity
          call PlMatCal (keycut, 
     &                   prop(1), time, dtime,
     &                   nTens, nlkword(1,1), dsdePl(1,1),etaEl(1),
     &                   wt1(1), wt2(1), Zee(1,1),
     &                   pleq, dpleq,
     &                   sigy, moist, A_MoistOffst, plstate)
      else
c --- this can only happen if cropt.eq.0
          dperr(1) = dble ( elemId )
          dperr(2) = dble ( matId )
          dperr(3) = dble ( cropt )
          call erhandler('CreepPlastDriver',5000,3,
     &     'No creep model specified for element %I' //
     &     ', material %I, model %I', dperr(1), ' ')
      endif
c
      if (plstate.eq.0) then
          goto 500
      elseif (plstate.lt.0) then
          goto 999
      endif

c --- update stresses, and plastic strains
      call v_move(etaEl(1),stress(1),nTens)
      if (keyplL .ne. 0 ) then
c --- this three lines below correct the plastic increment in case of softening behaviour
c     by calculating the elastic strain directly from the stresses (UPD: fixed elsewhere)
c          call v_move (eEl(1),wt2(1),nTens)
c          call MatOrthoInv (prop(1),stress(1),nTens,eEl(1))
c          call v1mv2_v1 (wt2(1),eEl(1),nTens)

          call v1mv2_v1 (eEl(1),wt2(1),nTens)
          call v1pv2_v1 (ePl(1),wt2(1),nTens)
      endif
      
      
      if (keycrL .ne. 0 ) then
          call v1mv2_v1(eEl(1),wt3(1),nTens)
          call v1pv2_v1(eCr(1),wt3(1),nTens)
      endif

      if (keyplL .gt. 0) then
          C(1) = ZERO
          do i = 1, nTens
             C(1) = C(1) + ( stress0(i) + stress(i) ) * wt2(i)
          enddo
          eqpl(1)  = - pleq
          eqpl(2)  = plWork_t + HALF * C(1)
      endif

c --- calculate the creep work
      if (keycrL .gt. 0) then
          C(1) = ZERO
          do i = 1, nTens
             C(1) = C(1) + ( stress0(i) + stress(i) ) * wt3(i)
          enddo
          eqcr(1)  = creq
          eqcr(2)  = crWork_t + HALF * C(1)
          matmisc(2) = sigy
      endif

      goto 600
  500 continue

c --- update stress in case of elastic/unloading
      call v_move (sigElp(1), stress(1), nTens)

  600 continue

c --- calculate the elastic work
      C(1) = ZERO
      do i = 1, nTens
         C(1) = C(1) + stress(i) * eEl(i)
      enddo
      sedEl = sedEl + HALF * C(1)

      if (keycrL .gt. 0) sedCr = sedCr + eqcr(2)
      if (keyplL .gt. 0) sedPl = sedPl + eqpl(2)

c --- check maximum plastic strain and creep ratio
      C(1) = ZERO
      call v_zero (Flstress(1),6)
      Flstress(1) = eEl(1)
      Flstress(4) = eEl(2)
      Flstress(5) = eEl(4)
      C(2) = EqvStrain (Flstress(1),HALF)
      if (C(2) .gt. ZERO .and. dpcr .gt. ZERO) then
          C(1) = dpcr / C(2)
      endif
c --- check maximum plastic strain and incremental  plastic strain
      sdv(1) = max(sdv(1), dpleq)
      sdv(2) = max(sdv(2), abs(pleq))
      sdv(3) = max(sdv(3), abs(creq))

 999  continue
c
      return
      end
