      subroutine PlMatCal (keycut,
     &                     prop, time, dtime,
     &                     nTens, nlkword, dsdePl, etaEl,
     &                     wt1, wt2, Zee,
     &                     pleq, dpleq, 
     &                     sigy, moist, MoistOffst, plstate)
c --- primary function: return mapping for creep with isotropic hardening for beam elements
c
c --- Reference
c           Based on Elastoplasticity and viscoplasticity
c           Simo, J.C. and Hughes, T.J.R., 1998
c
c  input arguments:
c      nTens    (int,sc,in)                -  nDirect + nShear
c      time     (dp,sc,in)                 -  current time
c      dtime    (dp,sc,in)                 -  current time increment
c      moist    (dp,sc,in)                 -  moisture at the matl point
c      propcr   (dp ,ar(*),in)             -  creep properties array
c      sigElp   (dp,ar(nTens),in)          -  stresses from elastic predictor
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
      double precision HALF,     ONE,      TWO,     SMALL,
     &                 ZERO,     TWOTHIRD, THREE,
     &                 SQTWOTHIRD, ONEDM02, SQONEHALF
      parameter       (ZERO       = 0.d0,
     &                 HALF       = 0.5d0,
     &                 ONE        = 1.d0,
     &                 TWO        = 2.d0,
     &                 THREE      = 3.d0,
     &                 SMALL      = 1.d-08,
     &                 ONEDM02    = 1.d-02, 
     &                 TWOTHIRD   = 2.0d0/3.0d0,
     &                 SQONEHALF  = 1.224744871391588900d0,      
     &                 SQTWOTHIRD = 0.816496580927726030d0,
     &                 NEWTON     = 15 ,
     &                 NDF        = 2)
c

c
      integer          plstate, elemId, nTens, matId, keycut
      integer*8        nlkword(W_nMatCol,W_nMatRow)

      double precision prop(*), etaEl(nTens), dtime, time,
     &                 wt1(nTens), wt2(nTens), wt4(nTens),
     &                 pleq, dpleq, sigyp, sigy, phi2, phi,
     &                 moist, dplpl, MoistOffst,
     &                 ds0dep, dsdePl(nTens,nTens)


      integer          kNewton, NewtonCont
      integer*8        keyplL

      double precision
     &                 pleq_t, twoG, sqTiny, 
     &                 dfdlm, lamb_pl, Fconsist, Fdef,
     &                 ResidR, dlamb_pl, qEl, Ex, G,
     &                 C(10), NC(10), two3rdE, Zee(nTens,nTens)
c
      double precision dperr(3), errNewton(NEWTON)
      logical          NewtonOK

      integer          option, nprop, plopt
      double precision proppl(48)

      keyplL = nlkword(1,1)
      plopt  = nlkword(2,1)
      
      plstate = 0
      ds0dep  = ZERO         
      sigyp   = sigy
      sqTiny  = sqrt(TINY)          
      lamb_pl = ZERO         
      dfdlm   = ZERO         
      NewtonCont = 0         

      pleq_t   = pleq
      dpleq    = ZERO
      Ex = prop(1)
      G  = prop(7)

      two3rdE = TWOTHIRD * Ex
      twoG    = Ex / (ONE + G)
c
      call v_move(etaEl(1), wt1(1), nTens)

      phi2 = TWOTHIRD * wt1(1) * wt1(1) + TWO * 
     &       (wt1(2) * wt1(2) +  wt1(3) * wt1(3))
c
      if (phi2 .le. ZERO) goto 999
      plstate = 1
      phi  = sqrt(phi2)
      qEl  = SQONEHALF * phi

c --- perform plasticity analysis
c
      if (plopt .eq. 3) then
          call GetMLIN (keyplL, moist, pleq, sigyp, ds0dep)
      else
          call GetTabParam(keyplL,moist,nprop,proppl(1),option)
          call YieldFun(nprop,prop(1),proppl(1),
     &                  qEl,pleq,dpleq,sigyp,ds0dep)
      endif
      
      if (sigyp .le. ZERO) goto 999

      ResidR = phi/sigyp - SQTWOTHIRD         
      if (ResidR.lt.-SMALL .or. abs(ResidR).lt.SMALL) then
          goto 999
      endif         

      plstate  = 2
      C(1) = TWOTHIRD * wt1(1) * wt1(1)
      C(2) = TWO * (wt1(2) * wt1(2) + wt1(3) * wt1(3))         
      C(3) = -C(1) * two3rdE
      C(4) = -C(2) * twoG          
      Fconsist = phi - SQTWOTHIRD * sigyp
      dfdlm = (C(3) + C(4)) / phi         
      Fdef = dfdlm - TWOTHIRD * ds0dep * phi
      dlamb_pl = -Fconsist / Fdef         
      lamb_pl = dlamb_pl         
      lamb_pl = max(lamb_pl, sqTiny)         
   10 continue                  
c         
      do kNewton = 1,NEWTON         
         NC(2) = ONE + two3rdE * lamb_pl         
         NC(1) = ONE + twoG * lamb_pl         
         NC(3) = ONE / NC(2)         
         NC(4) = NC(3) * NC(3)         
         NC(5) = ONE / NC(1)         
         NC(6) = NC(5) * NC(5)         
         phi2  = C(1) * NC(4) + C(2) * NC(6)         
         
         phi   = sqrt(phi2)         
         qEl   = SQONEHALF * phi
         dpleq = lamb_pl * SQTWOTHIRD * phi         
         pleq  = pleq_t + dpleq

         if (plopt .eq. 3) then
             call GetMLIN (keyplL, moist, pleq, sigyp, ds0dep)
         else
             call YieldFun(nprop,prop(1),proppl(1),
     &                     qEl,pleq,dpleq,sigyp,ds0dep)
         endif         
         
         Fconsist = phi - SQTWOTHIRD * sigyp         
         dfdlm    =  ( C(3) * NC(4) * NC(3)         
     &               + C(4) * NC(6) * NC(5)) / phi         
         dplpl    = SQTWOTHIRD * (lamb_pl * dfdlm + phi)         
         Fdef     = dfdlm - SQTWOTHIRD * ds0dep * dplpl          
         ResidR   = Fconsist/sigyp         
         dlamb_pl = -Fconsist / Fdef         
         lamb_pl  = lamb_pl + dlamb_pl           
         lamb_pl  = max(lamb_pl, sqTiny)         
         errNewton(kNewton) = abs(ResidR)

         if ((abs(ResidR).lt.SMALL .or. ResidR.lt.-SMALL)
     &        .and. ((abs(dlamb_pl) .lt. ONEDM02 * lamb_pl)
     &        .or. (dlamb_pl*lamb_pl .le. TINY))) goto 100             

      enddo      

c --- check whether it is worth to do Newton-raphson again         
      NewtonOK =         
     &   NewtonCont         .lt.4  .and.         
     & ((errNewton(NEWTON  ).lt.ONE.and.         
     &   errNewton(NEWTON  ).le.errNewton(NEWTON-1).and.         
     &   errNewton(NEWTON-1).le.errNewton(NEWTON-2)).or.         
     &  (errNewton(NEWTON  ).le.errNewton(NEWTON-1).and.         
     &   errNewton(NEWTON-1).le.errNewton(NEWTON-2).and.         
     &   errNewton(NEWTON-2).le.errNewton(NEWTON-3)))         
      NewtonCont = NewtonCont + 1
      
      if (NewtonOK) goto 10
      if (NewtonCont.lt.2.and.errNewton(NEWTON).le.ONEDM02) goto 10

      plstate  = - 1

c --- quit if uncovergence in plasticity algorithm
      dperr(1) = dble(plopt)
      call erhandler('PlMatCal',5000,2,
     &    'Integration algorithm did not converge for '//
     &    'Plastic model %I',dperr(1),' ')
      keycut = 1
      goto 999

 100  continue

c --- update stresses
      call v_zero(Zee(1,1),9)
      Zee(1,1) = THREE * Ex / (THREE + TWO * Ex * lamb_pl)
      Zee(2,2) = HALF * twoG / (lamb_pl * twoG + ONE)
      Zee(3,3) = Zee(2,2)

      wt1(1) = Zee(1,1) * wt1(1) / Ex
      wt1(2) = Zee(2,2) * wt1(2) * TWO / twoG
      wt1(3) = Zee(3,3) * wt1(3) * TWO / twoG
      call v_move(wt1(1), etaEl(1), nTens)

c --- form the consistent elasto-plastic moduli
      wt4(1) = Zee(1,1) * TWOTHIRD * wt1(1)
      wt4(2) = Zee(2,2) * TWO * wt1(2)
      wt4(3) = Zee(3,3) * TWO * wt1(3)

      C(1) = TWOTHIRD * wt4(1) * wt1(1) + 
     &         TWO * (wt4(2) * wt1(2) + wt4(3) * wt1(3))
      C(2) = TWOTHIRD * ds0dep * phi2 / 
     &       ( ONE - TWOTHIRD * ds0dep * lamb_pl )
      C(1) = ONE / ( C(1) + C(2) )

      dsdePl(1,1) = Zee(1,1) - (C(1) * wt4(1) * wt4(1))  
      dsdePl(2,2) = Zee(2,2) - (C(1) * wt4(2) * wt4(2))   
      dsdePl(3,3) = Zee(3,3) - (C(1) * wt4(3) * wt4(3))
      dsdePl(1,2) = ZERO 
      dsdePl(1,3) = ZERO
      dsdePl(2,3) = ZERO
      dsdePl(2,1) = dsdePl(1,2)
      dsdePl(3,1) = dsdePl(1,3)
      dsdePl(3,2) = dsdePl(2,3)

c --- update incremental plastic strain
      wt2(1) = TWOTHIRD * lamb_pl *  wt1(1)
      wt2(2) = TWO * lamb_pl * wt1(2)
      wt2(3) = TWO * lamb_pl * wt1(3)

 999  continue
      return
      end