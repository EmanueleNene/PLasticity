      subroutine ChabocheSV (matId, elemId, keycut,
     &                       ncomp, nStatev, nProp,
     &                       stress, ustatev, dsdePl, sedEl, sedPl, 
     &                       Strain, dStrain, epsPl, prop)
c*************************************************************************
c     *** primary function ***
c     von Mises material with nonlinear kinematic hardening by Chabouche
c
c*************************************************************************
c
c     input arguments
c     ===============
c      matId     (int,sc,i)               material #
c      elemId    (int,sc,i)               element #
c      nStatev   (int,sc,l)               Number of state variables
c      nProp     (int,sc,l)               Number of material ocnstants
c
c      Strain   (dp,ar(ncomp),i)          Strain at beginning of time increment
c      dStrain  (dp,ar(ncomp),i)          Strain increment
c      prop     (dp,ar(nprop),i)          Material constants defined by TB,USER
c
c     input output arguments              
c     ======================             
c      stress   (dp,ar(nTesn),io)         stress
c
c      ustatev   (dp,ar(nStatev),io)       statev
c           ustatev(1) - ustatev(ncomp)    - plastic strain vector
c           ustatev(4)                     - equivalent plastic strain
c           ustatev(5) - ustatev(4+ncomp)  - backstress vector
c           ustatev(5+ncomp)               - elastic state (1 - el; 0 - pl)   
c           ustatev(nStatev)               - von-Mises stress
c
c      sedEl    (dp,sc,io)                elastic work
c      sedPl    (dp,sc,io)                plastic work
c
c     output arguments
c     ================
c      keycut   (int,sc,io)               loading bisect/cut control
c                                         0 - no bisect/cut
c                                         1 - bisect/cut 
c                                         (factor will be determined by ANSYS solution control)
c      dsdePl   (dp,ar(ncomp,ncomp),io)   material jacobian matrix
c
c      stresss and strains, plastic strain vectors
c          11, 13, 12                for 3d beam
c      
c
c      material jacobian matrix
c        3d beam plasticity
c           dsdePl    |  1111   1113   1112 |
c           dsdePl    |  1311   1313   1312 |
c           dsdePl    |  1211   1213   1212 |
c      
c     elastic constants:
c         young    = prop(1)
c         posn     = prop(2)
c         sigy0    = prop(3)
c         R0       = prop(4)
c         Rinf     = prop(5)
c         b        = prop(6)
c         lame2 ==> mu
c*************************************************************************
#include "impcom.inc"
c
      integer          matId, elemId,
     &                 keycut,
     &                 ncomp,nStatev,nProp
      double precision
     &                 sedEl,   sedPl
      double precision
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 elstate, akin(ncomp),
     &                 bs(ncomp), dsigydga
c
c***************** User defined part *************************************
c
c --- parameters
c
      integer          NEWTON, mcomp
      double precision HALF, ONE, TWO, SMALL, SQTWOTHIRD,ONETHIRD,
     &                 ZERO, TWOTHIRD, ONEDM02, ONEDM05, sqTiny,
     &                 TOLER
      parameter       (ZERO       = 0.d0,
     &                 HALF       = 0.5d0,
     &                 ONE        = 1.d0,
     &                 TWO        = 2.d0,
     &                 SMALL      = 1.d-08,
     &                 sqTiny     = 1.d-10,
     &                 ONEDM02    = 1.d-02,
     &                 ONEDM05    = 1.d-05,
     &                 TWOTHIRD   = 2.0d0/3.0d0,
     &                 ONETHIRD   = 1.0d0/3.0d0,
     &                 SQTWOTHIRD = 0.816496580927726030d0,
     &                 NEWTON     = 30,
     &                 mcomp      = 3,
     &                 TOLER      = 1.0d-5
     &                 
     &                       )
c
c --- local variables
c
c      sigElp   (dp,ar(3  ),l)            trial stress
c      dsdeEl   (dp,ar(3,3),l)            elastic moduli
c      pleq_t   (dp,sc     ,l)            equivalent plastic strain at beginnig of time increment
c      pleq     (dp,sc     ,l)            equivalent plastic strain at end of time increment
c      dpleq    (dp,sc     ,l)            incremental equivalent plastic strain
c      gamma    (dp,sc     ,l)            variable for solving incremental equivalent plastic strain
c      dgamma   (dp,sc     ,l)            correction of gamma
c      sigy_t   (dp,sc     ,l)            yield stress at beginnig of time increments
c      sigy     (dp,sc     ,l)            yield stress at end of time increment
c      young    (dp,sc     ,l)            Young's modulus
c      posn     (dp,sc     ,l)            Poiss's ratio
c      sigy0    (dp,sc     ,l)            initial yield stress
c      dsigdep  (dp,sc     ,l)            plastic slop
c      twoG     (dp,sc     ,l)            two time of shear moduli
c      threeG   (dp,sc     ,1)            three times shear moduli 
c      funcf    (dp,sc     ,l)            nonlinear function to be solved for gamma
c      dFdep    (dp,sc     ,l)            derivative of nonlinear function over gamma
c      akin     (dp,ar(3  ),1)            kinematic hardening backstress (alpha)
c      dsdevdga (dp,ar(3  ),1)            deviatoric stress derivative w.r.t. gamma
c      dsigydga (dp, sc    ,1)            derivative of sigy w.r.t. gamma

c --- temperary variables for solution purpose
c      i, j
c      c1, c2, c3, fratio, den, abseta
c      wk1(3), wk2(3), wk3(3), wk4(3) vector working arrays
c
                    
      double precision sigElp(mcomp), dsdeEl(mcomp,mcomp), 
     &                 wk1(3), wk2(3), wk3(3), wk4(3), dsdevdga(mcomp),
     &                 sdev(mcomp),
     &                 M3(mcomp,mcomp),normal(mcomp),eta(mcomp),
     &                 eta_red(mcomp)

      double precision R0, Rinf, b, g1, CC1, kog, CC2

      integer          i, j, k
      double precision pleq_t,  sigy_t , sigy,
     &                 dpleq,   pleq,    twoG,
     &                 young, posn,    sigy0,   dsigdep, 
     &                 gamma, dgamma,  dfdga,   dplga,   fratio,
     &                 funcFb,funcFb2, funcf,   dFdep, threeG,
     &                 c1, c5, gammaK,
     &                 lame2,den,VM,young0, abseta,
     &                 dFdep2, dF, error

      double precision kn(3)
      data kn/ONE,      TWO, TWO/

c*************************************************************************
      keycut   = 0
c *** equivalent plastic strain at beginning of time step
      pleq_t   = ustatev(4)
      pleq     = pleq_t
      elstate = ustatev(5+ncomp)
c     VM - von Mises equivalent stressfrom prev state
      VM = ustatev(nStatev)       
c *** get Young's modulus and Poisson's ratio, initial yield stress and slope of stress-strain
      young0   = prop(1)
      posn     = prop(2)
      sigy0    = prop(3)
      R0       = prop(4)
      Rinf     = prop(5)
      b        = prop(6)
      g1       = prop(7)
      CC1      = prop(8)
      CC2      = prop(9)

c *** calculate plastic slope
      dsigdep  = 2*Rinf*pleq+ R0
c *** updated Young's modulus      
      young    = young0 + CC2 * pleq
      twoG     = young / (ONE+posn)
      threeG   = twoG/TWOTHIRD
      lame2    = twoG / TWO
      
c *** calculate elastic stiffness matrix
c ===============

      call v_zero(dsdeEl(1,1), ncomp * ncomp)
      c1 = twoG * HALF
      dsdeEl (1,1) = young
      dsdeEl (2,2) = c1
      dsdeEl (3,3) = c1
c  wk3 = [E, G, G]       
      do i = 1, ncomp
         wk3(i) = dsdeEl(i,i)
      enddo
    
c *** calculate prescribed strain (by increment dStrain) eq.[2.1]
c *** [wk1] = [Strain]+[dStrain]-[ustatev(1:3)]
      call v_move(Strain(1), wk1(1), ncomp)
      call v1pv2_v1(wk1(1), dStrain(1), ncomp)
      call v1mv2_v1(wk1(1), ustatev(1), ncomp)
      
c *** get state variables from previous step (backstress)
      call v_zero(akin(1),ncomp)
      call v_move(ustatev(5), akin(1), ncomp)
c *** bs - backstress at step n - needed for BS update      
      call v_move(ustatev(5), bs(1), ncomp)

c     copy elastic moduli dsdeEl to material Jacobian matrix dsdePl
      call v_move(dsdeEl(1,1), dsdePl(1,1), ncomp * ncomp)

c *** calculate the trial stress [sigElp] = [dsdeEl].[wk1] eq.[2.6]
      do i=1,ncomp
         sigElp(i) = wk3(i) * wk1(i)
      enddo

c *** calculate deviatoric stress of the trial stress eq.[2.5]
      call v_zero(sdev(1),ncomp)
      call v_move(sigElp(1), sdev(1), ncomp)
      do i=1,ncomp
         sdev(i) = sdev(i)-ONETHIRD*sigElp(i)
      enddo
c     sdev(1) = sdev(1)- ONETHIRD*sigElp(1)
      
c *** calculate working stress [eta] = [sdev] - [akin]   eq.[3.9]    
      do i=1,ncomp
         eta(i) = sdev(i)-akin(i)
      enddo
      
c *** calculate von Mises stress - trial eq.[2.6]
      funcFb2 = ZERO
      do i = 1, ncomp
         funcFb2 = funcFb2 + kn(i)*eta(i) * eta(i)
      enddo
c *** funcFb - von Mises equivalent stress VM = sqrt(3*(eta11^2+2*(eta12^2+eta13^2))/2)      
      funcFb = sqrt(funcFb2/TWOTHIRD)/SQTWOTHIRD
      if (funcFb .le. ZERO) goto 990

c *** compute current yield stress
c      sigy    = sigy0 + Rinf*exp(b*pleq)
      sigy = sigy0 + Rinf*pleq*pleq + R0*pleq
c --- eq.[2.7]       
      fratio = funcFb/sigy - ONE

c      IF (elstate .eq. ZERO) THEN
c     if the previous step was elastic check for yielding  eq.[2.8]    
          if (fratio .LE. -SMALL) then 
              goto 500
          else
              goto 50
          endif    
c      ENDIF                  
      
   50 continue
      elstate = ONE
                  
c     PLASTIC STEP
      ustatev(5+ncomp) = ONE
      sigy_t = sigy

      gammaK = ZERO
      dplga = ZERO
      dfdga = ZERO
      dpleq = ZERO
      pleq = pleq_t 

c *** Local New-Raphson procedure for solving the gamma [Appendix B]
c     (plastic multiplier) and thus, the incremental 
c     equivalent plastic strain (dlpeq)
      do k = 1,NEWTON
         gamma = gammaK
         funcFb2 = ZERO
         dfdga = ZERO
c ***    n+1 yield stress eq.[6.1]
c         sigy  = sigy0 + Rinf*exp(b*(pleq+gamma))
         sigy = sigy0 + Rinf*(pleq+gamma)*(pleq+gamma)+R0*(pleq+gamma)
c ----  eq.[2.10]                                        
         call v_zero(eta(1),ncomp)                                    
         kog = ONE/(ONE + g1*gamma)
         if (kog .gt. ONE) then
             call erhandler('ChabocheSV',5000,2,
     &       'PROBLEMS IN DEVIATORIC STRESS',ZERO,' ')             
             keycut = 1
             goto 990      
         endif
c ----  eq.[6.2]
         do i = 1, ncomp
            den = sigy + (threeG + kog * CC1) * gamma
            eta(i) = sigy*(sdev(i) - kog*akin(i))/den
            eta_red(i) = sdev(i) - kog*akin(i)
         enddo

c        calculate von Mises equiv stress(reduced) eq.[6.4]
c        von Mises equivalent stress VM = sqrt(3*(eta_red11^2+2*(eta_red12^2+eta_red13^2))/2)      
         do i = 1, ncomp
            funcFb2 = funcFb2 + kn(i) * eta_red(i) * eta_red(i)
c           dsdevdga = d(eta_red)/(d(gamma)            
            dsdevdga(i) = kog*kog*g1*akin(i)
            dfdga = dfdga + kn(i) * eta_red(i) * dsdevdga(i)
         enddo
         funcFb = sqrt(funcFb2/TWOTHIRD)/SQTWOTHIRD 
         
c ***    derivative of funcFb w.r.t. gamma (dfdga*3/2) eq.[6.5]
         dfdga = (dfdga/TWOTHIRD)/SQTWOTHIRD 

c ***    derivative of sigy(n+1) w.r.t. gamma: dsigydga=d(sigy)/d(gamma)   eq.[6.9]
c         dsigydga = Rinf * b * exp(b*(pleq+gamma))
         dsigydga = Rinf*TWO*(pleq+gamma) + R0

c ***    RESIDUAL::::::::::::::::::::::::::::::::
c ***    calculate the residual eq.[6.7]
         funcf = funcFb - sigy - gamma*(threeG + kog*CC1)
c ***    derivative of residual function w.r.t. gamma dFep = d(funcf)/d(gamma) eq.[6.8]
         dFdep = dfdga/funcFb - dsigydga -(threeG - kog*CC1)
         dFdep2 = kog*kog*g1*CC1
         dF = dFdep + gamma*dFdep2
              
c *** NEWTON SOLUTION FOR GAMMA::::::::::::::::
c *** correction of gamma
         dgamma = -funcf / dF 
         gammaK = gamma   + dgamma
c *** check for negative gamma eq.[6.10]
         gammaK = max (gammaK, sqTiny)
         error = abs (gammaK - gamma)
         fratio = funcf/ sigy
        
c *** Check for convergence of local New-Raphson iteration
  
         if (((abs(fratio) .lt. ONEDM05 ) .AND.
     &        (abs(dgamma) .lt. ONEDM02*gammaK)) .OR.
     &       ((abs(fratio) .lt. ONEDM05 ) .AND.
     &        (dgamma      .le. sqTiny  ) .AND.
     &        ( gammaK     .le. sqTiny  )))  goto 100

      enddo

c *** check for UNconvergence myself
      if (error .ge. TOLER) then
          call erhandler('ChabocheSV',5000,2,
     &    'PLASTIC ALGORITH DID NOT CONVERGE',ZERO,' ')
          keycut = 1
      endif

      if (gammaK .lt. ZERO) then
          call erhandler('ChabocheSV',5000,2,
     &    'NEGATIVE PLASTIC MULTIPLIER',ZERO,' ')
          keycut = 1
      endif
c *** Uncovergence, set keycut to 1 for bisection/cutback the load increment
      keycut = 1
      goto 990
          
  100 continue
      gamma = gammaK
c *** calculate the incremental equivalent plastic strain eq.[3.7]
      dpleq = gamma                               
c *** update the total equivalent plastic strain eq.[3.8]
      pleq = pleq_t + dpleq

c *** UPDATE ALL THE VALUES
c *** updated Young's and shear modulus  eq.[3.9]     
      young = young0 + CC2 * pleq
      twoG = young / (ONE+posn)
c *** update elastic stiffness matrix      
      call v_zero(dsdeEl(1,1), ncomp * ncomp)
      c1 = twoG * HALF
      dsdeEl (1,1) = young
      dsdeEl (2,2) = c1
      dsdeEl (3,3) = c1
c *** copy elastic moduli dsdeEl to material Jacobian matrix dsdePl
      call v_move(dsdeEl(1,1), dsdePl(1,1), ncomp * ncomp)
c *** update yield stress (n+1) eq.[2.11]
c      sigy     = sigy0 + Rinf*exp(b*(pleq))
      sigy = sigy0 + Rinf*pleq*pleq + R0*pleq
c *** current plastic slope
      dsigdep  = 2*Rinf*pleq+ R0
      
c *** update eta (n+1)  eq.[3.10] 
c *** abseta = ||eta|| - absolute value of a vector   eq.[3.12]       
      kog = ONE/(ONE + g1*gamma)
      abseta  = 0
      do i = 1, ncomp
         den = sigy + (threeG + kog * CC1) * gamma
         eta(i) = sigy*(sdev(i) - kog*akin(i))/den
         abseta = abseta + eta(i)*eta(i)
      enddo
      abseta = sqrt(abseta)

c *** calculate normal to the yield surface [normalized vector n(i) = eta(i)/||eta||] eq.[3.11]
      do i = 1, ncomp
          normal(i) = eta(i)/(abseta/SQTWOTHIRD)
      enddo 
            
c *** calculate incremental plastic strain wk2=epsPl eq.[3.13]
      do j = 1, ncomp
         wk2(j) = gamma * normal(j)/SQTWOTHIRD 
      enddo
c *** update plastic strains eq.[3.14]
      call v1pv2_v1 (epsPl(1),wk2(1),ncomp)
c *** calculate stresses   eq.[3.15]  
      call v_zero(wk4(1),ncomp)
      den = 0
      do i = 1, ncomp
         do j = 1, ncomp
            den = den + dsdeEl(i,j)*wk2(j)
         enddo
         wk4(i) = sigElp(i) - den
         den = 0
      enddo

c *** save stresses wk4
      call v_move(wk4(1), stress(1), ncomp) 
      
c *** update back stresses eq.[3.15]
      kog = ONE/(ONE + g1*gamma)
      do i = 1, ncomp
         akin(i) = kog*(bs(i)+TWOTHIRD*CC1*wk2(i))
      enddo  
      
c *** Update state variables
      ustatev(ncomp+1) = pleq
      
      do i=1,ncomp
         ustatev(i) = epsPl(i)
      enddo

      do i=5,4+ncomp
         ustatev(i) = akin(i-4)
      enddo
      
c *** update plastic work
      sedPl = sedPl + HALF * (sigy_t+sigy) * dpleq
      
      call v_zero(M3(1,1), ncomp * ncomp)
c *** c5 - slope eq.[4.1]
      c5 = TWO*TWO*lame2*lame2/(TWO*lame2+sigy)
      do i = 1, ncomp
         do j = 1, ncomp
            M3(i,j) = c5*normal(i)*normal(j)
         enddo
      enddo   
      
c *** Calculate the consistent Jacobian eq.[4.1]
      do i=1,ncomp
         do j=1,ncomp
            dsdePl(i,j) = dsdePl(i,j)- M3(i,j)
         enddo
      enddo

      goto 600

  500 continue
   
c *** Update stress in case of elastic/unloading eq.[2.8]
      call v_move (sigElp(1),stress(1),ncomp)
      elstate = ZERO
  600 continue
c *** elastic strain energy
      sedEl = ZERO
      do i = 1, ncomp
         sedEl = sedEl + stress(i)*(Strain(i)+dStrain(i)-epsPl(i))
      enddo
      sedEl = sedEl * HALF
      ustatev(nStatev) = funcFb
      ustatev(5+ncomp) = elstate

  990 continue
          
      return
      end