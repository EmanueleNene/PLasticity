      subroutine ChabocheCalc (keycut, matId, nlkword, prop,
     &                         time, dtime, 
     &                         nTens, plstate,
     &                         numxi, sigy0, ckh, gkh,
     &                         xback0, xback, xbacki, dsdePl, 
     &                         etaEl, wt1, wt2, wt3, wt4, 
     &                         dii, pleq, dpleq, sigyp, 
     &                         moist, plWork)
c *** primary function: energy calculations for beam element
c                       using return mapping algorithm for 
c                       rate-independent Chaboche plasticity model
c
#include "impcom.inc"
#include "constants.inc"
c
      integer          NEWTON, NUMF, maxhn
      double precision HALF, ONE, TWO, SMALL, ONEHALF,
     &                 ZERO, ONEDM02, ONEDM05, TWOTHIRD,
     &                 SQTWOTHIRD, SQONEHALF
      parameter       (ZERO = 0.0d0,
     &                 HALF = 0.5d0,
     &                 ONE = 1.d0,
     &                 TWO = 2.d0,
     &                 SMALL = 1.d-08,
     &                 ONEDM02 = 1.d-02,
     &                 ONEDM05 = 1.d-05,
     &                 ONEHALF = 1.5d0,
     &                 TWOTHIRD = 2.0d0/3.0d0,
     &                 SQONEHALF = 1.224744871391588900d0,
     &                 SQTWOTHIRD = 0.816496580927726030d0,
     &                 NEWTON = 15,
     &                 NUMF = 1)

      integer          keycut, npropISO, optISO,
     &                 plstate, matId, nTens, numxi
      double precision dii(nTens), dsdePl(nTens,nTens),
     &                 prop(*), propISO(1000),
     &                 xback0(nTens), xback(nTens), 
     &                 xbacki(nTens,numxi),etaEl(nTens), dtime, time,
     &                 wt1(nTens), wt2(nTens), wt3(nTens), wt4(nTens),
     &                 plWork, pleq, dpleq, sigyp, sigy, moist

      integer          i, j, kNewton, NewtonCont
      double precision con1, con2, con3, con4,
     &                 twoG , oneOv3G, tolSigy, sqTiny, 
     &                 alamda, adlamda, dfdga,
     &                 funcFb, funcFb1, funcFb2, dplga,
     &                 dgamma, gamma, fratio, pleq_t, 
     &                 threeG, ds0dep, sigy0, sigy_t, sigy1,
     &                 ci, gi
      double precision gkh(numxi), ckh(numxi)
c
      double precision dperr(3), f(NUMF), fs(NUMF,NUMF),
     &                 errNewton(NEWTON)
      double precision pv(3)
      data             pv /TWOTHIRD, TWO, TWO/
      logical          NewtonOK
      integer*8        nlkword(W_nMatCol,W_nMatRow), keyplL, keycpL
c
      plstate = 0
      kNewton = 0

      ds0dep = ZERO
      call v_zero (wt2(1),nTens)

      twoG = prop(1) / (ONE + prop(7))
      threeG = ONEHALF * twoG
      oneOv3G = ONE / threeG
c
      funcFb2 = ZERO
      do i = 1,nTens
         funcFb2 = funcFb2 + pv(i) * etaEl(i) * etaEl(i)
      enddo
      funcFb = sqrt(funcFb2)

c
      if (funcFb.le.ZERO) goto 900

      tolSigy = SMALL * oneOv3G
      sqTiny = sqrt (TINY) / funcFb
      sigy = SQONEHALF * funcFb
      sigy_t = sigyp
c
      gamma = ZERO
      sigy1 = sigy0

      keyplL = nlkword(1,1)
      keycpL = nlkword(1,3)
c --- isotropic hardening plasticity
      if (keycpL.ne.keyplL) then
          call GetTabParam (keyplL, moist, npropISO, propISO(1), optISO)
          call GetIsoYield (propISO(1),optISO,pleq,sigy1,ds0dep,sigy0)
      endif
          
      sigyp = sigy1

      if (sigyp.lt.ZERO) goto 900
      fratio = sigy/sigyp - ONE
      if (fratio.lt.-SMALL) then
          pleq = - pleq
          goto 900
      endif
      if (sigy_t.le.ZERO) sigy_t = sigyp

      do i = 1,nTens
         dii(i) = dsdePl(i,i)
      enddo

      do i = 1,nTens
         dii(i) = dii(i) * pv(i)
      enddo

      plstate = 1
      pleq_t = pleq 

      NewtonCont = 0
      dpleq = gamma * SQTWOTHIRD * funcFb
      pleq = pleq_t + dpleq
      dplga = ZERO
      funcFb1 = ONE / funcFb

   10 continue
c
      do kNewton = 1,NEWTON

         alamda = ZERO
         adlamda = ZERO
         call v_zero (wt2(1), nTens)
         call v_zero (wt3(1), nTens)
         do i = 1,numxi
            ci = ckh(i)
            gi = gkh(i)
            if (abs(gi).gt.ZERO) then
                con1 = exp(- gi * dpleq)
                con2 = ONE - con1
                alamda = alamda  + con2 * ci / gi
                adlamda = adlamda + ci * con1
                do j = 1,nTens
                   wt2(j) = wt2(j) + con2 * xbacki(j,i)
                enddo
                con2 = con1 * gkh(i)
                do j = 1,nTens
                   wt3(j) = wt3(j) + con2 * xbacki(j,i)
                enddo
             else
                alamda = alamda  + ci * dpleq
                adlamda = adlamda + ci
             endif
         enddo

         con1 = ONE / sigyp
         alamda = con1 * alamda
         adlamda = con1 * adlamda
         adlamda = adlamda - alamda * con1 * ds0dep

         call v1pv2_v1 (wt2(1), etaEl(1), nTens)

         do j = 1,nTens
            con1 = ONE + alamda + gamma * dii(j)
            con1 = ONE / con1
            wt1(j) = con1
            con2 = wt2(j) * con1
            wt2(j) = con2
            con3 = con2 * con1
            wt4(j) = wt3(j) * con1 - con3 * adlamda
         enddo
         funcFb2 = ZERO
         do j = 1,nTens
            con1 = wt2(j)
            funcFb2 = funcFb2 + pv(j) * con1 * con1
         enddo
         
         funcFb = sqrt(funcFb2)
         funcFb1 = ONE / funcFb

         con1 = ZERO
         do j = 1,nTens
            con1 = con1 + pv(j) * wt2(j) * wt4(j)
         enddo
         con2 = ZERO
         do j = 1,nTens
            con2 = con2 + pv(j) * wt2(j) * wt1(j) * dii(j) * wt2(j)
         enddo

         dplga = TWOTHIRD * sigyp / 
     &          (ONE - TWOTHIRD * gamma * ds0dep)
         dfdga = funcFb1 * (con1 * dplga - con2)
         f(1) = funcFb - SQTWOTHIRD * sigyp
         fratio = abs(f(1))/sigyp
         fs(1,1) = dfdga - SQTWOTHIRD * ds0dep * dplga

         dgamma = -f(1) / fs(1,1)
         gamma = gamma + dgamma  
         gamma = max (gamma, sqTiny)
         dpleq = TWOTHIRD * gamma * sigyp
         pleq = pleq_t + dpleq
         if (keycpL.ne.keyplL) then
             call GetIsoYield(propISO(1),optISO,pleq,sigy1,ds0dep,sigy0)
         endif    
         sigyp = sigy1
         errNewton(kNewton) = abs(fratio)
c        

         if (((abs(fratio) .lt. ONEDM05 ) .and.
     &        (abs(dgamma) .lt. ONEDM02*gamma)) .or.
     &       ((abs(fratio) .lt. ONEDM05 ) .and.
     &        (abs(gamma ) .le. tolSigy )))  goto 100

      enddo

c --- check whether it is worth to do Newton-raphson again
      NewtonOK =
     &   NewtonCont         .lt.4  .and.
     & ((errNewton(NEWTON  ).lt.ONE.and.
     &   errNewton(NEWTON  ).lt.errNewton(NEWTON-1).and.
     &   errNewton(NEWTON-1).lt.errNewton(NEWTON-2)).or.
     &  (errNewton(NEWTON  ).lt.errNewton(NEWTON-1).and.
     &   errNewton(NEWTON-1).lt.errNewton(NEWTON-2).and.
     &   errNewton(NEWTON-2).lt.errNewton(NEWTON-3)))
      NewtonCont = NewtonCont + 1

      if (NewtonOK) goto 10
      if (NewtonCont.lt.2 .and. errNewton(NEWTON).le.ONEDM02) goto 10

c --- uncovergence in plasticity algorithm
      plstate = - 1
      dperr(1) = dble(matId )
      call erhandler('ChabocheCalc',5000,2,
     &    'Plasticity integration failed to converge for '//
     &    'material %I, Cycling plasticty',dperr(1),' ')

      keycut = 1
      goto 900

 100  continue

c --- update back stress
      sigy = funcFb * SQONEHALF
      call v_move (wt2(1), etaEl(1), nTens)
      call v_zero (xback(1), nTens)
      do i = 1,numxi
         ci = ckh(i)
         gi = gkh(i)
         if (abs(gi).gt.ZERO) then
             con1 = exp(- gi * dpleq)
             con2 = SQTWOTHIRD * ci / gi * funcFb1
             con3 = ONE - con1
             do j = 1,nTens
                con4 = (con2 * etaEl(j) - xbacki(j,i)) * con3
                xbacki(j,i) = xbacki(j,i) + con4
             enddo
         else
             con1 = SQTWOTHIRD * ci * funcFb1 * dpleq
             do j = 1,nTens
                xbacki(j,i) = xbacki(j,i) + con1 * etaEl(j)
             enddo
         endif
         call v1pv2_v1 (xback(1), xbacki(1,i), nTens)
      enddo
c --- wt2(i) stores incremental plastic strain and shouldn't be reused
      do j = 1, nTens
         wt2(j) = gamma * pv(j) * etaEl(j)
      enddo
      con1 = ZERO
      do j = 1, nTens
         con1 = con1 + (xback0(j) + xback(j)) * wt2(j)
      enddo
      plWork = plWork + HALF * con1 + HALF * (sigy_t + sigyp) * dpleq

      adlamda  = adlamda * dplga
      call cv1_v1 (wt3(1), nTens, dplga)
c
      con1 = (ONE + alamda)
      con2 = TWOTHIRD * ds0dep
      con4 = con2 * funcFb2 / (ONE - con2 * gamma)
      con4 = con4 * con1 + adlamda * funcFb2
      do j = 1,nTens
         con4 = con4 - pv(j) * etaEl(j) * wt3(j)
      enddo
      con1 = ONE / con1
      con2 = con1 * gamma
      con3 = ONE - con2 * adlamda
      call v_move(etaEl(1), wt1(1), nTens)
      do j = 1,nTens
         wt4(j) = wt1(j) * con3 + con2 * wt3(j)
      enddo
      do j = 1,nTens
         con3 = ONE / (ONE + con2 * dii(j))
         dii(j) = dii(j) * con3 / pv(j)
      enddo
      do j = 1,nTens
         wt1(j) = wt1(j) * pv(j)
      enddo
      do j = 1,nTens
         wt4(j) = wt4(j) * pv(j)
      enddo
      do j = 1,nTens
         con4 = con4 + wt1(j) * wt4(j) * dii(j)
      enddo
      do j = 1,nTens
         wt1(j) = wt1(j) * dii(j)
      enddo
      do j = 1,nTens
         wt4(j) = wt4(j) * dii(j)
      enddo

      con4 = ONE / con4
      do i = 1,nTens
         dsdePl(i,i) = dii(i)
      enddo
      do i = 1,nTens
         do j = 1,nTens
            dsdePl(i,j) = dsdePl(i,j) - con4 * wt4(i) * wt1(j)
         enddo
      enddo
c --- symmetrilize the consistent elastoplastic tangent moduli
      do i = 1,nTens - 1
         do j = i+1,nTens
            con1 = HALF*(dsdePl(i,j)+dsdePl(j,i))
            dsdePl(i,j) = con1
            dsdePl(j,i) = con1
         enddo
      enddo

 900  continue
c
      return
      end