      subroutine ContDebond (nvmax,cinem,qa,qb,presc0,taurc0,tausc0,
     &                       dfunc,rigp0,rigt0a,rigt0s,idbmode,stiff,
     &                       ONED)
c --- primary function: calculate stiffness for debonded contact elements
c
c  input arguments:
c       nvmax   (int,sc,in)         - number of state variables
c       cinem   (dp,ar(39),in)      - kinematic parameters and material data
c                                     used for contact and debonding
c       qa      (dp,ar(*),in)       - stresses and state variables at initial
c                                     configuration
c       qb      (dp,ar(*),in)       - stresses and state variables at updated
c                                     configuration
c       presc0  (dp,sc,in)          - normal contact pressure w/o damage
c       taurc0  (dp,sc,in)          - tangent contact stress in r-dir w/o damage
c       tausc0  (dp,sc,in)          - tangent contact stress in s-dir w/o damage
c       dfunc   (dp,sc,in)          - damage function value at current substep
c       rigp0   (dp,sc,in)          - initial normal penalty stiffness
c       rigt0   (dp,sc,in)          - initial tangent penalty stiffness
c       idbmode (int,sc,inout)      - failure mode for debonding
c
c  output arguments:
c       stiff   (dp,ar(4,4),out)    - stiffness matrix
c
#include "impcom.inc"
#include "solu_info.inc"
c
      integer nvmax,idbmode
c
      double precision cinem(41),qa(nvmax),qb(nvmax),stiff(4,4),
     &  presc0,taurc0,tausc0,dfunc,rigp0,rigt0a,rigt0s

      integer icflg
c
      double precision ZERO,ONE,HALF,ONED,ASMALL,dtest2,
     &                 tnmax,ttmax,gnmax,gtmax,vrpar,tauac0,
     &                 dnew,dold,dvnew,sfratio,fact1,fact2,vfac1
c
      parameter (ZERO=0.0d0,ONE=1.0d0,HALF=0.5d0,ASMALL=1.d-8)
      
#if defined(DEBUG)
      call SubBeg('ContDebond')
#endif

c --- initialize
c
      sfratio = ZERO
      vfac1 = ZERO
      tauac0 = sqrt(taurc0*taurc0 + tausc0*tausc0)! Combined tangent traction

c --- cohesive zone model data for debonding
c
      tnmax = cinem(36)   ! Max normal traction
      gnmax = cinem(37)
      ttmax = cinem(38)   ! Max tangent traction
      gtmax = cinem(39)
      vrpar = cinem(40)
      icflg = nint(cinem(41))

c --- retrieve damage variables
c
      dnew = qb(nvmax-9)
      dvnew = qb(nvmax-11)
      dold = qa(nvmax-9)

c --- assemble stiffness
c
      if (dvnew.gt.ZERO .and. dvnew.lt.ONED) then !damage
          if (presc0.lt.ZERO) then
              stiff(3,3) = rigp0*(ONE-dvnew)
              stiff(1,1) = rigt0a*(ONE-dvnew)
              stiff(2,2) = rigt0s*(ONE-dvnew)
          else
              stiff(3,3) = rigp0
              if ((idbmode.eq.1) .or.
     &            (idbmode.eq.3 .and. icflg.eq.0)) then
                   stiff(1,1) = rigt0a
                   stiff(2,2) = rigt0s
              else
                   stiff(1,1) = rigt0a*(ONE-dvnew)
                   stiff(2,2) = rigt0s*(ONE-dvnew)
              endif
          endif
          stiff(1,3) = ZERO
          stiff(2,3) = ZERO
          stiff(3,1) = ZERO
          stiff(3,2) = ZERO
          stiff(1,2) = ZERO
          stiff(2,1) = ZERO
          dtest2 = abs((dnew/dold) - ONE)
          if (dnew.gt.dold .and. dtest2.gt.ASMALL) then ! new damage
              vfac1 = A_TimeInc/(A_TimeInc+vrpar)
              if (idbmode.eq.1) then
                  sfratio = ONE - tnmax/(rigp0*gnmax)
                  fact1 = sfratio*abs(presc0)
                  fact2 = tnmax*vfac1
                  stiff(3,3) = stiff(3,3) - rigp0*fact2/fact1
              elseif (idbmode.eq.2) then
                  sfratio = ONE - ttmax/(rigt0a*gtmax)
                  fact1 = sfratio*tauac0
                  fact2 = ttmax*vfac1
                  stiff(1,1) = stiff(1,1) - rigt0a*fact2/fact1
                  sfratio = ONE - ttmax/(rigt0s*gtmax)
                  fact1 = sfratio*tauac0
                  fact2 = ttmax*vfac1
                  stiff(2,2) = stiff(2,2) - rigt0s*fact2/fact1
              elseif (idbmode.eq.3) then
                  if (presc0.lt.ZERO) then
                      sfratio = ONE - tnmax/(rigp0*gnmax)
                      fact1 = sfratio*(dfunc+ONE)**3
                      fact2 = vfac1*presc0*presc0/(tnmax*tnmax)
                      stiff(3,3) = stiff(3,3) - rigp0*fact2/fact1
                      fact2 = vfac1*taurc0*taurc0/(ttmax*ttmax)
                      stiff(1,1) = stiff(1,1) - rigt0a*fact2/fact1
                      fact2 = vfac1*tausc0*tausc0/(ttmax*ttmax)
                      stiff(2,2) = stiff(2,2) - rigt0s*fact2/fact1
                      fact2 = vfac1*presc0*taurc0/(ttmax*ttmax)
                      stiff(3,1) = stiff(3,1) - rigt0a*fact2/fact1
                      fact2 = vfac1*presc0*tausc0/(ttmax*ttmax)
                      stiff(3,2) = stiff(3,2) - rigt0s*fact2/fact1
                      fact2 = vfac1*presc0*taurc0/(tnmax*tnmax)
                      stiff(1,3) = stiff(1,3) - rigp0*fact2/fact1
                      fact2 = vfac1*presc0*tausc0/(tnmax*tnmax)
                      stiff(2,3) = stiff(2,3) - rigp0*fact2/fact1
c --- symmetrize
                      fact2 = HALF*(stiff(1,3)+stiff(3,1))
                      stiff(1,3) = fact2
                      stiff(3,1) = fact2
                      fact2 = HALF*(stiff(2,3)+stiff(3,2))
                      stiff(2,3) = fact2
                      stiff(3,2) = fact2
                 else
                      sfratio = ONE - tnmax/(rigp0*gnmax)
                      fact1 = sfratio*(dfunc+ONE)**3
                      fact2 = vfac1*taurc0*taurc0/(ttmax*ttmax)
                      stiff(1,1) = stiff(1,1) - rigt0a*fact2/fact1
                      fact2 = vfac1*tausc0*tausc0/(ttmax*ttmax)
                      stiff(2,2) = stiff(2,2) - rigt0s*fact2/fact1
                 endif
              endif
          endif
      elseif (dvnew.ge.ONED) then ! failure
          stiff(1,1) = ZERO
          stiff(1,2) = ZERO
          stiff(2,1) = ZERO
          stiff(2,2) = ZERO
          stiff(1,3) = ZERO
          stiff(3,1) = ZERO
          stiff(3,3) = ZERO
      elseif (dvnew.eq.ZERO) then ! no damage
          stiff(1,1) = rigt0a
          stiff(1,2) = ZERO
          stiff(2,1) = ZERO
          stiff(2,2) = rigt0s
          stiff(1,3) = ZERO
          stiff(3,1) = ZERO
          stiff(3,3) = rigp0
      endif

#if defined(DEBUG)
      call SubEnd('ContDebond')
#endif
      
      return
      end
