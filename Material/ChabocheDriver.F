      subroutine ChabocheDriver (keycut, time, dtime, moist, dmoist,
     &                           nTens, elemId, matId, prop,
     &                           nxiMax, propcp0, propcp,
     &                           eEl, ePl, eqpl, 
     &                           xback, xbacki, matmisc,
     &                           wt1, wt2, wt3, wt4,
     &                           etaEl, xback0, stress0, wt21,
     &                           sigElp, stress, dsdePl, 
     &                           nlkword, sedEl, sedPl, sdv)
c --- primary function: cyclic plasticity driver
c
c  input arguments
c      keycut   (int,sc,in)                - bisec/cut flag
c      moist    (dp,sc,in)                 - moisture at the mat point
c      nTens    (int,sc,in)                - nDirect + nShear
c      sigElp   (dp,ar(nTens),in)          - elastic predictor
c      nlkword  (int,ar(2,i),in)           - nonlinear material index
c
c  output arguments
c      eEl      (dp,ar(nTens),inout)       - eEl_t + eIncTot/ eEl_tDt
c      ePl      (dp,ar(nTens),inout)       - ePl_t          / ePl_tDt
c      stress   (dp,ar(nTesn),inout)       - stress(elprd)  / stress_tDt
c      dsdePl   (dp,ar(nTens,nTens),inout) - el stiff       / el-pl stiff
c      eqpl     (dp,ar(*),inout)           - eqpl_t         / eqpl_tDt
c                                          - plWork_t       / plWork_tDt
c      xback    (dp,ar(nTens),inout)       - xback_t        / xback_tDt
c
c  working arrays and local variables
c      wt1      (dp,ar(*),work)            - work tensor - Sel
c      wt2      (dp,ar(*),work)            - work tensor - flow dir
c      wt3      (dp,ar(*),work)            - work tensor - inc pl strain
c      wt21     (dp,ar(nTens,*),work)      - work tensor - 4th order I'
c
c      sigy_t   (dp,sc,l)                  - yield stress at 't'
c      gmu      (dp,sc,l)                  - shear moduli
c      eEldev   (dp,ar(6),l)               - deviatoric elastic strain at t+Dt
c      etaEl    (dp,ar(6),l)               - deviatoric back-stress    at t+Dt
c
#include "impcom.inc"
#include "constants.inc"
c
      double precision ZERO, HALF, TOLER
      
      parameter       (ZERO = 0.0d0,
     &                 HALF = 0.5d0,
     &                 TOLER = 1.0d-8)
      
      integer*8        nlkword(W_nMatCol,W_nMatRow)
      integer          keycut,
     &                 nTens, 
     &                 elemId, matId,
     &                 numxi, nxiMax
      double precision dsdePl(nTens,nTens),
     &                 eEl(nTens), ePl(nTens), eqpl(2), 
     &                 sigElp(nTens), etaEl(nTens), 
     &                 xbacki(nTens,nxiMax),
     &                 wt1(nTens), wt2(nTens), wt3(nTens), wt4(nTens),
     &                 wt21(nTens,nTens), xback0(nTens), stress0(nTens),
     &                 stress(nTens), sedEl, sedPl,
     &                 xback(nTens), sdv(*), matmisc(*)
      double precision time, dtime, moist, dmoist
      double precision propcp0(*), propcp(*)

      integer          i, plstate
      integer*8        keycpL
      double precision plWork, sigy, pleq, dpleq, 
     &                 prop(13), ci, gi, ci0, gi0


      keycpL = nlkword(1,3)

      plstate = 0
      dpleq = ZERO 
      sedEl = ZERO
      sedPl = ZERO
      call v_zero (xback(1), nTens)

      call GetChabProp (keycpL, moist, numxi, 
     &                  propcp(1), propcp(2),
     &                  propcp(2+nxiMax), propcp(2+2*nxiMax))
c --- in case of moisture change, we need to calculate the contribution
c     due to rate of change of Ci. note, xbacki is the increment of back stresses
      if (abs(dmoist).gt.TOLER) then
          call GetChabProp (keycpL, moist-dmoist, numxi, 
     &                      propcp0(1), propcp0(2), 
     &                      propcp0(2+nxiMax), propcp0(2+2*nxiMax))
      else
          call v_move (propcp(1),propcp0(1),1+nxiMax*3)
      endif
      
      do i = 1,numxi
         ci = propcp(1+i)
         gi = propcp(1+nxiMax+i)
         ci0 = propcp0(1+i)
         gi0 = propcp0(1+nxiMax+i)
         if (abs(gi0).gt.ZERO) then
c --- this accounts for both Ci and Gi
             ci = (ci * gi0)/(gi * ci0)
             call cv1_v1 (xbacki(1,i), nTens, ci)
         elseif (abs(ci0).gt.ZERO) then
c --- if there is no Gi specified
             ci = ci / ci0
             call cv1_v1 (xbacki(1,i), nTens, ci)
         endif
         call v1pv2_v1 (xback(1), xbacki(1,i), nTens)
      enddo      

      call v_move (xback(1), xback0(1), nTens)
      call v_move (stress(1), stress0(1), nTens)
      pleq = abs(eqpl(1))
      plWork = eqpl(2)
      sigy = matmisc(2)

      call v1mv2_v3 (sigElp(1), xback0(1), etaEl(1), nTens)

      call ChabocheCalc (keycut, matId, nlkword(1,1), prop(1),
     &                   time, dtime,
     &                   nTens, plstate,
     &                   numxi, propcp(1), propcp(2),
     &                   propcp(2+nxiMax),
     &                   xback0(1), xback(1), xbacki(1,1), 
     &                   dsdePl(1,1), 
     &                   etaEl(1), wt1(1), wt2(1), wt3(1),
     &                   wt4(1), wt21(1,1), 
     &                   pleq, dpleq, sigy, moist, plWork)

      if (plstate.eq.0) then
c --- elastic / unloading
          eqpl(1) = - eqpl(1)
          goto 500
      elseif (plstate.lt.0) then
c --- plasticity algorihm unconverge
          goto 990
      endif

c --- check maximum plastic strain and incremental plastic strain
      sdv(1) = max(sdv(1), dpleq)
      sdv(2) = max(sdv(2), abs(pleq))

c --- update back-stresses, plastic strains and stresses

      call v1pv2_v3 (etaEl(1), xback(1), stress(1), nTens)
      call v1pv2_v1 (ePl(1), wt2(1), nTens)
      call v1mv2_v1 (eEl(1), wt2(1), nTens)

      eqpl(1) = - pleq
      eqpl(2) = plWork

      goto 600
  500 continue

c --- update stress in case of elastic/unloading
      call v_move (sigElp(1), stress(1), nTens)

  600 continue
c --- elastic strain energy
      do i = 1,nTens
         sedEl = sedEl + stress(i) * eEl(i)
      enddo
      sedEl = HALF * sedEl
      sedPl = eqpl(2)
      matmisc(2) = sigy
      matmisc(1) = sedEl

 990  continue

      return
      end