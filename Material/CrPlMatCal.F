      subroutine CrPlMatCal (keycut, elemId, matId,
     &                       prop, propcr, time, dtime,
     &                       nTens, nlkword, dsdePl, etaEl,
     &                       wt1, wt2, wt3, Zee,
     &                       pleq, creq, dpleq, dpcr,
     &                       sigy, moist, MoistOffst, plstate)
c --- primary function: return mapping for creep with isotropic hardening for beam elements
c
c --- Reference
c           Based on Elastoplasticity and viscoplasticity
c           Simo, J.C. and Hughes, T.J.R., 1998
c
c  input arguments:
c      nTens    (int,sc,in)                -  nDirect + nShear
c      time                                -  current time
c      dtime                               -  current time increment
c      moist    (dp,sc,in)                 -  moisture at the matl point
c      propcr   (dp ,ar(*),in)             -  creep properties array
c      sigElp   (dp,ar(nTens),in)          -  sresees from elastic predictor
c
c  output arguments:
c      sigy     (dp,sc,out)                -  yield stress/pressure at 't' / 't+dt'
c      plstate                             -  plastic status 0 - not; 1 - yes
c      eEl      (dp,ar(nTens),inout)       -  eEl_t + eIncTot/ eEl_tDt
c      ePl      (dp,ar(nTens),inout)       -  ePl_t / ePl_tDt
c      dsdePl   (dp,ar(nTens,nTens),inout) -  el stiff / el-pl stiff
c      eqpl     (dp,ar(*),inout)           -  eqpl_t / eqpl_tDt
c                                          -  plWork_t / plWork_tDt
c      wt2      (dp,ar(*),inout)           -  incremental plastic strain
c      wt3      (dp,ar(*),inout)           -  incremental creep strain
c
#include "impcom.inc"
#include "constants.inc"
c
      integer          NEWTON, NDF
      double precision HALF,     ONE,      TWO,     SMALL,    ONEHALF,
     &                 ZERO,     TWOTHIRD, THREE,   ONEDM01,  ONEDM03,
     &                 SQTWOTHIRD, SQONEHALF, FOURONINE, FOUR
      parameter       (ZERO       = 0.d0,
     &                 HALF       = 0.5d0,
     &                 ONE        = 1.d0,
     &                 TWO        = 2.d0,
     &                 THREE      = 3.d0,
     &                 FOUR       = 4.d0,
     &                 SMALL      = 1.d-08,
     &                 ONEDM01    = 5.d-01,
     &                 ONEDM03    = 1.d-03,
     &                 ONEHALF    = 1.5d0,
     &                 TWOTHIRD   = 2.0d0/3.0d0,
     &                 SQONEHALF  = 1.224744871391588999d0,
     &                 SQTWOTHIRD = 0.816496580927726030d0,
     &                 FOURONINE  = 4.0d0/9.0d0,
     &                 NEWTON     = 15 ,
     &                 NDF        = 2)
c

c
      integer     plstate, elemId, nTens, matId, keycut
      integer*8   nlkword(W_nMatCol,W_nMatRow)

      double precision propcr(*), prop(*), etaEl(nTens), dtime, time,
     &                 wt1(nTens), wt2(nTens), wt3(nTens), wt4(nTens),
     &                 pleq, creq, dpleq, dpcr, sigyp, sigy,
     &                 moist, dcr, dplpl, dplcr, dcrpl, MoistOffst,
     &                 ds0dep, dsdePl(nTens,nTens)



      integer          kNewton, NewtonCont, cropt
      integer*8        keyplL

      double precision
     &                 pleq_t, creq_t, con1, con2, con3, con4,
     &                 twoG, sqTiny, Fcnsst0, funcFb, Fcnsst2, det,
     &                 dfdlm, gammacr, lamb_pl, dpcrT, lambda,
     &                 dcreq, TResult, FTResult, dlamb_cr, dlamb_pl,
     &                 C(10), NC(10), phi2, phi,
     &                 two3rdE, fplast, Zee(nTens,nTens)
c
      double precision dcrds(4), dperr(3), errNewton(NEWTON),
     &                 F(2), fs(2,2), epsref
      logical          NewtonOK

      integer          option, nprop
      double precision proppl(48), sigy0, Et

      lamb_pl  = ZERO
      gammacr  = ZERO
      dlamb_cr = ZERO
      ds0dep   = ZERO
      fplast   = ZERO
      sigyp    = ZERO
      dcreq    = ZERO


      keyplL = nlkword(1,1)
      cropt  = nlkword(2,2)

      plstate = 0
c
      two3rdE = TWOTHIRD * prop(1)
      twoG    = prop(1) / (ONE + prop(7))
      call v_move(etaEl(1), wt1(1), nTens)
c
      phi2 = TWOTHIRD * wt1(1) * wt1(1) + TWO *
     &      (wt1(2) * wt1(2) +  wt1(3) * wt1(3) )
      phi  = sqrt ( phi2 )
      if (phi .le. ZERO ) goto 999
      if (sigy .le. ZERO) sigy = phi
c
      call v_zero(dcrds(1),4)

      creq_t = creq
      creq = creq_t + dcreq
      sqTiny = sqrt(TINY)
      plstate = 1
      funcFb = sigy

      NewtonCont = 0
c
      call CreepLaws (cropt, propcr(1), time, dtime,
     &                moist, MoistOffst, creq_t,
     &                dpcr, dcrds(1), creq, sigy)

      dlamb_cr = dpcr / ( SQTWOTHIRD * sigy )
      gammacr  = gammacr + dlamb_cr

c --- perform pure creep analysis
c
      con1 = TWOTHIRD * wt1(1) * wt1(1)
      con2 = two3rdE
      con3 = TWO * (wt1(2) * wt1(2) + wt1(3) * wt1(3))
      con4 = twoG
      Fcnsst0 = funcFb
	epsref = sigy/twoG

  10  continue
c
      do kNewton = 1, NEWTON
         call CreepFun ( cropt, propcr(1), time, dtime,
     &                   moist, MoistOffst,
     &                   creq_t, dpcr, dcrds(1), dpcrT,
     &                   creq, sigy,
     &                   con1, con2, con3, con4,
     &                   gammacr, funcFb, NDF, f(1), fs(1,1), epsref)

         FTResult  = (funcFb - Fcnsst0) / funcFb
         Fcnsst0  = funcFb
         TResult   = max ( abs(FTResult), abs(F(1)) )
c
         dlamb_cr = -F(1) / fs(1,1)
         gammacr  = gammacr + dlamb_cr
         gammacr  = max( gammacr, sqTiny )
         errNewton(kNewton) = abs(TResult)

         if ((abs(TResult) .lt. ONEDM03) .and.
     &      ((abs(dlamb_cr) .lt. ONEDM03*gammacr) .or.
     &      ((abs(dlamb_cr) .lt. SMALL) .and.
     &       (abs(gammacr) .lt. SMALL))))  goto 100
      enddo

c --- check whether it is worth to do Newton-raphson again
      NewtonOK =
     &   NewtonCont         .lt.4  .and.
     & ((errNewton(NEWTON  ).lt.ONE.and.
     &   errNewton(NEWTON  ).le.ONEDM01*errNewton(NEWTON-1).and.
     &   errNewton(NEWTON-1).le.ONEDM01*errNewton(NEWTON-2)).or.
     &  (errNewton(NEWTON  ).le.ONEDM01*errNewton(NEWTON-1).and.
     &   errNewton(NEWTON-1).le.ONEDM01*errNewton(NEWTON-2).and.
     &   errNewton(NEWTON  ).le.ONEDM03*errNewton(NEWTON-8)))
      NewtonCont = NewtonCont + 1
      if (NewtonOK) goto 10
c
      plstate = - 1
c --- uncovergence in creep algorithm
      dperr(1) = dble(cropt)
      call erhandler('CrPlMatCal',5000,2,
     &    'Integration algorithm did not converge for '//
     &    'Material model %I',dperr(1),' ')
      keycut = 1
      goto 999
 100  continue
      lambda = gammacr

      if (keyplL.eq.0) goto 200

c --- perform combined creep and plasticity analysis
c
      call GetTabParam(keyplL,moist,nprop,proppl(1),option)
      sigy0 = abs(proppl(1))
      Et = abs(proppl(2))
      ds0dep = prop(1) * Et / (prop(1) - Et)
      sigyp = sigy0 + pleq * prop(1) * Et / (prop(1) - Et)

      if (sigy .lt. sigyp) then
          pleq = - pleq
          goto 200
      endif

      if (ds0dep.le.ZERO) ds0dep = 1.0E-20
      plstate  = 2
      pleq_t = pleq
      fplast = ONE
      NewtonCont = 0

 110  continue
c
      con1 = TWOTHIRD * wt1(1) * wt1(1)
      con2 = two3rdE
      con3 = TWO * (wt1(2) * wt1(2) + wt1(3) * wt1(3))
      con4 = twoG

      do kNewton = 1, NEWTON

         lambda = lamb_pl + gammacr
         C(1) = ONE + two3rdE * lambda
         C(2) = ONE + twoG * lambda
         NC(1) = ONE / C(1)
         NC(2) = NC(1) * NC(1)
         NC(3) = ONE / C(2)
         NC(4) = NC(3) * NC(3)
         Fcnsst2 = con1 * NC(2) + con3 * NC(4)
         funcFb = sqrt ( Fcnsst2 )
c
         dcreq = gammacr * SQTWOTHIRD * funcFb
         dpleq = lamb_pl * SQTWOTHIRD * funcFb
         creq = creq_t + dcreq
         pleq = pleq_t + dpleq
         sigy = funcFb * SQONEHALF

         call CreepLaws (cropt, propcr(1), time, dtime,
     &                   moist, MoistOffst, creq_t,
     &                   dpcr, dcrds(1), creq, sigy)
c
         ds0dep = prop(1) * Et / (prop(1) - Et)
         sigyp = sigy0 + pleq * ds0dep
c --- consistent slopes

         F(1) = dpcr - SQTWOTHIRD * gammacr * funcFb
         F(2) = funcFb - SQTWOTHIRD * sigyp

         dfdlm = - (con1 * con2 * NC(2) * NC(1)
     &           +  con3 * con4 * NC(4) * NC(3) ) / funcFb
         dcr = SQTWOTHIRD * (gammacr * dfdlm + funcFb)
         dplpl = SQTWOTHIRD * (lamb_pl * dfdlm + funcFb)
         dplcr = SQTWOTHIRD * lamb_pl * dfdlm
         dcrpl = SQTWOTHIRD * gammacr * dfdlm
         fs(1,1) = SQONEHALF * dcrds(1) * dfdlm +
     &               dcrds(3) * dcr - dcr
         fs(1,2) = SQONEHALF * dcrds(1) * dfdlm +
     &               dcrds(3) * dcrpl - dcrpl
         fs(2,1) = dfdlm - SQTWOTHIRD * ds0dep * dplcr
         fs(2,2) = dfdlm - SQTWOTHIRD * ds0dep * dplpl

         det = fs(1,1) * fs(2,2) - fs(2,1) * fs(1,2)
         det = - ONE / det
         dlamb_cr = det * (fs(2,2) * F(1) - fs(1,2) * F(2))
         dlamb_pl = det * (fs(1,1) * F(2) - fs(2,1) * F(1))
         gammacr = gammacr + dlamb_cr
         lamb_pl = lamb_pl + dlamb_pl
         gammacr = max (gammacr, sqTiny)

         TResult = F(1)
         if (dpcr .gt. SMALL ) TResult = abs( F(1) / dpcr )
         if (lamb_pl .gt. sqTiny ) then
             TResult = max(abs(ONE - SQTWOTHIRD * sigyp/funcFb),TResult)
         else
             F(2) = ZERO
             fs(1,2) = ZERO
             fs(2,1) = ZERO
             lamb_pl = sqTiny
         endif
         errNewton(kNewton) = abs(TResult)
         FTResult = (funcFb - Fcnsst0) / funcFb
         Fcnsst0 = funcFb
         TResult = max ( abs(FTResult), abs(TResult) )
c
         if (((abs(TResult)   .lt. ONEDM03         ) .and.
     &        (abs(dlamb_pl) .lt. ONEDM03*lamb_pl ) .and.
     &        (abs(dlamb_cr) .lt. ONEDM03*gammacr)) .or.
     &       ((abs(TResult)   .lt. ONEDM03 ) .and.
     &        (abs(dlamb_pl) .le. SMALL ) .and.
     &        (abs(dlamb_cr) .le. SMALL ) .and.
c     &        (abs(lamb_pl ) .le. SMALL ) .and.
     &        (abs(gammacr ) .le. SMALL )))  goto 200
      enddo

c --- check whether it is worth to do Newton-raphson again
      NewtonOK =
     &   NewtonCont         .lt.4  .and.
     & ((errNewton(NEWTON  ).lt.ONE.and.
     &   errNewton(NEWTON  ).le.ONEDM01*errNewton(NEWTON-1).and.
     &   errNewton(NEWTON-1).le.ONEDM01*errNewton(NEWTON-2)).or.
     &  (errNewton(NEWTON  ).le.ONEDM01*errNewton(NEWTON-1).and.
     &   errNewton(NEWTON-1).le.ONEDM01*errNewton(NEWTON-2).and.
     &   errNewton(NEWTON  ).le.ONEDM03*errNewton(NEWTON-8)))
      NewtonCont = NewtonCont + 1
      if (NewtonOK) goto 110

      plstate  = - 1
c --- quit if uncovergence in plasticity algorithm
      dperr(1) = dble(cropt)
      call erhandler('CrPlMatCal',5000,2,
     &    'Integration algorithm did not converge for '//
     &    'Creep model %I',dperr(1),' ')
      keycut = 1
      goto 999

 200  continue

c --- update stresses
c --- during the update, lambda = lamb_pl + gammacr
      call v_zero(Zee(1,1),9)
      Zee(1,1) = THREE * prop(1) / (THREE + TWO * prop(1) * lambda)
      Zee(2,2) = HALF * twoG / (lambda * twoG + ONE)
      Zee(3,3) = Zee(2,2)

      wt1(1) = Zee(1,1) * wt1(1) / prop(1)
      wt1(2) = Zee(2,2) * wt1(2) * TWO / twoG
      wt1(3) = Zee(3,3) * wt1(3) * TWO / twoG
      
      call v_move(wt1(1), etaEl(1), nTens)      

c --- form the consistent elasto-creep-plastic moduli

      Fcnsst2 = funcFb * funcFb
      con1 = - (ONEHALF * dcrds(1) + (dcrds(3)-ONE) * gammacr) /
     &         (Fcnsst2 * (dcrds(3)-ONE) )
      if (ds0dep .le. 1.d-18) ds0dep = 1.d-18
      con2 = fplast * (ONEHALF - lamb_pl * ds0dep) /
     &                (Fcnsst2 * ds0dep)
      con3 = ONE / ( con1 + con2 )
      con2 = FOURONINE * Zee(1,1) * wt1(1) * wt1(1) +
     &           Zee(2,2) * FOUR * (wt1(2) * wt1(2) + wt1(3) * wt1(3))

      con1 = ONE / ( con2 + con3 )

      wt4(1) = Zee(1,1) * TWOTHIRD * wt1(1)
      wt4(2) = Zee(2,2) * TWO * wt1(2)
      wt4(3) = Zee(3,3) * TWO * wt1(3)

      dsdePl(1,1) = Zee(1,1) - (con1 * wt4(1) * wt4(1))
      dsdePl(2,2) = Zee(2,2) - (con1 * wt4(2) * wt4(2))
      dsdePl(3,3) = Zee(3,3) - (con1 * wt4(3) * wt4(3))
      dsdePl(1,2) = ZERO
      dsdePl(1,3) = ZERO
      dsdePl(2,3) = ZERO
      dsdePl(2,1) = dsdePl(1,2)
      dsdePl(3,1) = dsdePl(1,3)
      dsdePl(3,2) = dsdePl(2,3)

c --- update incremental creep strain
      wt3(1) = TWOTHIRD * gammacr * wt1(1)
      wt3(2) = TWO * gammacr * wt1(2)
      wt3(3) = TWO * gammacr * wt1(3)

c --- update incremental plastic strain
      if (keyplL.gt.0) then
c --- update incremental plastic strain
          wt2(1) = TWOTHIRD * lamb_pl *  wt1(1)
          wt2(2) = TWO * lamb_pl * wt1(2)
          wt2(3) = TWO * lamb_pl * wt1(3)
      else
          call v_zero(wt2(1), nTens)
      endif
c
 900  continue

 999  continue
      return
      end
