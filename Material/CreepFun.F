      subroutine CreepFun (cropt, propcr, time, dtime,
     &                     moist, MoistOffst,
     &                     creq_t, dpcr, dcrds, dpcrT,
     &                     creq, sigy, con1, con2, con3, con4,
     &                     gammacr, funcFb, NDF, f, fs, epsref)
c --- primary function: calculate value of creep functions and its derivitives
c
#include "impcom.inc"
c
      double precision ONE, SQTWOTHIRD, SQONEHALF
      parameter       (ONE        = 1.0d0,
     &                 SQONEHALF  = 1.224744871391588900d0,
     &                 SQTWOTHIRD = 0.816496580927726030d0)

      integer          NDF, cropt
      double precision f(NDF),fs(NDF,NDF), propcr(*),
     &                 time, dtime, moist, MoistOffst,
     &                 creq_t, dpcr, dcrds(*),
     &                 funcFb, gammacr, creq, sigy,
     &                 con1, con2, con3, con4
      double precision con2p1, ocon2p1,
     &                 ocon2p2, con4p1, ocon4p1, ocon4p2,
     &                 dfdga, dacrcr ,dcrdga, dpcrT, epsref,
     &                 SMALLCR
      logical          isCreqSmall

      con2p1 = ONE + con2 * gammacr
      con4p1 = ONE + con4 * gammacr
      ocon2p1 = ONE / con2p1
      ocon2p2 = ocon2p1 * ocon2p1
      ocon4p1 = ONE / con4p1
      ocon4p2 = ocon4p1 * ocon4p1
      funcFb = sqrt(con1 * ocon2p2 + con3 * ocon4p2)
      dpcrT = SQTWOTHIRD * gammacr * funcFb
      creq = creq_t + dpcrT
      sigy = funcFb * SQONEHALF
      call CreepLaws (cropt, propcr(1), time, dtime,
     &                moist, MoistOffst, creq_t,
     &                dpcr, dcrds(1), creq, sigy)

      dfdga  =  - (con1 *con2 * ocon2p2 * ocon2p1
     &          +  con3 * con4 * ocon4p2 * ocon4p1)/funcFb
      dacrcr = SQTWOTHIRD * (gammacr * dfdga + funcFb)
      dcrdga = SQONEHALF * dcrds(1) * dfdga + dcrds(3)*dacrcr

      SMALLCR = 1.0d-3 * epsref
      isCreqSmall = (abs(dpcr).eq.0) .or. (abs(dpcrT).lt.SMALLCR)

      if (isCreqSmall) then
          f(1) = dpcr - dpcrT
          fs(1,1) = SQONEHALF * dcrds(1) * dfdga +
     &              (dcrds(3)-ONE) * dacrcr
      else
          f(1) = ONE - dpcrT/dpcr
          fs(1,1) = - dacrcr/dpcr + dpcrT*dcrdga/(dpcr*dpcr)
      endif


      return
      end
