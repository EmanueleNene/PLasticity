      subroutine CrMatCal (keycut,
     &                     prop, propcr, time, dtime,
     &                     nTens, nlkword, dsdePl, etaEl,
     &                     wt1, wt3, Zee,
     &                     creq, dpcr,
     &                     sigy, moist, MoistOffst, Crstate)
c --- primary function: return mapping for pure creep for beam elements
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
c      wt3      (dp,ar(*),inout)           -  incremental creep strain
c
#include "impcom.inc"
#include "constants.inc"
c
      integer          NEWTON, NDF
      double precision HALF, THIRD, ONE, TWO, SMALL, ONEHALF,
     &                 ZERO, ONEDM02, ONEDM03, TWOTHIRD, THREE, 
     &                 SQTWOTHIRD, FOUR, FOURONINE,
     &                 SQONEHALF
      parameter       (ZERO       = 0.d0,
     &                 HALF       = 0.5d0,
     &                 THIRD      = 1.d0/3.d0,
     &                 ONE        = 1.d0,
     &                 TWO        = 2.d0,
     &                 THREE      = 3.d0,
     &                 FOUR       = 4.d0,
     &                 SMALL      = 1.d-10,
     &                 ONEDM02    = 1.d-03,
     &                 ONEDM03    = 1.d-04,
     &                 ONEHALF    = 1.5d0,
     &                 TWOTHIRD   = 2.0d0/3.0d0,
     &                 FOURONINE  = 4.0d0/9.0d0,
     &                 SQONEHALF  = 1.224744871391588999d0,
     &                 SQTWOTHIRD = 0.816496580927726030d0,
     &                 NEWTON     = 25,
     &                 NDF        = 1)

c --- argument list

      integer          Crstate, nTens
      integer          keycut 
      integer*8        nlkword(W_nMatCol,W_nMatRow)

      double precision propcr(*), dsdePl(nTens,nTens),
     &                 prop(*), etaEl(*), dtime, time, ZPeta(3),
     &                 wt1(*), wt3(*), Zee(3,*),
     &                 creq, dpcr, dpcrT,
     &                 sigy, moist, MoistOffst, test_n

      integer          kNewton, cropt
      double precision pres,   creq_t,
     &                 con1,   con2,   con3,   con4,   
     &                 twoG,   tolSigy, sqTiny, 
     &                 funcFb0, funcFb, funcFb2, 
     &                 dcreq, dgammacr, fratio, gammacr, 
     &                 phi2, phi, two3rdE, epsref

      double precision dcrds(4), dperr(3), f(1), fs(1,1)

      Crstate = 0
      cropt   = abs(nlkword(2,2))
      kNewton = 0
      call v_zero(dcrds(1),4)
      gammacr = ZERO

      two3rdE = TWOTHIRD * prop(1)
      twoG    = prop(1) / (ONE + prop(7))
      call v_move( etaEl(1), wt1(1), nTens )

      phi2 = TWOTHIRD * wt1(1) * wt1(1) + TWO * 
     &      (wt1(2) * wt1(2) +  wt1(3) * wt1(3) ) 
      phi = sqrt(phi2)

      if (phi .le. ZERO) goto 999
      if (sigy.le. ZERO) sigy = phi
      tolSigy = SMALL
      sqTiny  = sqrt(TINY)
      creq_t  = creq 
      dcreq = ZERO 
      pres  = ZERO 
      creq  = creq_t + dcreq
                       
      call CreepLaws (cropt, propcr(1), time, dtime,
     &                moist, MoistOffst, creq_t,
     &                dpcr, dcrds(1), creq, sigy)
c     
      dgammacr = dpcr / ( SQTWOTHIRD * sigy )
      gammacr  = gammacr + dgammacr

c --- perform pure creep analysis

      con1    = TWOTHIRD * wt1(1) * wt1(1)
      con2    = two3rdE
      con3    = TWO * ( wt1(2) * wt1(2) + wt1(3) * wt1(3) )
      con4    = twoG
      funcFb0 = sigy
	epsref = sigy/twoG

      do kNewton=1,NEWTON

         call CreepFun ( cropt, propcr(1), time, dtime,
     &                   moist, MoistOffst,
     &                   creq_t, dpcr, dcrds(1), dpcrT,
     &                   creq, sigy,
     &                   con1, con2, con3, con4,
     &                   gammacr, funcFb, NDF, f(1), fs(1,1), epsref)
c     
         fratio = (funcFb - funcFb0)/funcFb
         fratio = max(abs(fratio),abs(f(1)))
         funcFb0  = funcFb
c     
         if ((abs(fratio) .lt. ONEDM03 ) .and. 
     &      ((abs(dgammacr) .lt. ONEDM03*gammacr ) .or.   
     &      ((abs(dgammacr) .lt. SMALL ) .and.   
     &       (abs(gammacr ) .lt. SMALL))))  goto 100   
     
c --- new newton prediction
         dgammacr = -f(1) / fs(1,1)
         gammacr  = gammacr + dgammacr
         gammacr  = max(gammacr, sqTiny)

      enddo

c --- relax the convergence condition
      if (fratio .lt. tolSigy*10 .and. 
     &    abs(dgammacr) .lt. ONEDM02*gammacr) goto 100
      
      Crstate  = - 1
c --- uncovergence in creep algorithm is detected
      dperr(1) = abs(cropt)
      call erhandler('CrMatCal',5000,2,
     &    'Integration algorithm did not converge for '//
     &    'Material model %I',dperr(1),' ')
      keycut = 1
      goto 999
 100  continue

c --- update creep strains and stresses increments
      Crstate = 1
      call v_zero(Zee(1,1),9)
      Zee(1,1) = THREE * prop(1) / (THREE + TWO * prop(1) * gammacr)
      Zee(2,2) = HALF * twoG / (gammacr * twoG + ONE)
      Zee(3,3) = Zee(2,2)

      wt1(1) = Zee(1,1) * wt1(1) / prop(1)          
      wt1(2) = Zee(2,2) * wt1(2) * TWO / twoG 
      wt1(3) = Zee(3,3) * wt1(3) * TWO / twoG 
      call v_move(wt1(1), etaEl(1), nTens)
        
c --- update incremental creep strain
      wt3(1) = TWOTHIRD * gammacr * wt1(1)
      wt3(2) = TWO * gammacr * wt1(2)
      wt3(3) = TWO * gammacr * wt1(3)
      
      funcFb2 = funcFb * funcFb
      test_n = (ONEHALF * dcrds(1) + (dcrds(3)-ONE) * gammacr) 
      if (test_n .ne. ZERO) then
          con1 = -( funcFb2 * (dcrds(3)-ONE)) / test_n
      else
          con1 = -( funcFb2 * (dcrds(3)-ONE)) / TINY
      endif
      con2  = FOURONINE * Zee(1,1) * wt1(1) * wt1(1) +
     &        Zee(2,2) * FOUR * (wt1(2) * wt1(2) + wt1(3) * wt1(3))
      con1 = ONE / ( con1 + con2 )
      ZPeta(1) = Zee(1,1) * TWOTHIRD * wt1(1)
      ZPeta(2) = Zee(2,2) * TWO * wt1(2)
      ZPeta(3) = Zee(3,3) * TWO * wt1(3)
         
      dsdePl(1,1) = Zee(1,1) - (con1 * ZPeta(1) * ZPeta(1))  
      dsdePl(2,2) = Zee(2,2) - (con1 * ZPeta(2) * ZPeta(2))   
      dsdePl(3,3) = Zee(3,3) - (con1 * ZPeta(3) * ZPeta(3))
      dsdePl(1,2) = ZERO 
      dsdePl(1,3) = ZERO
      dsdePl(2,3) = ZERO
      dsdePl(2,1) = dsdePl( 1,2 )  
      dsdePl(3,1) = dsdePl( 1,3 )   
      dsdePl(3,2) = dsdePl( 2,3 )  
c     
 999  continue
      return
      end
