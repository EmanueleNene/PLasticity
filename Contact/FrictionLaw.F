      subroutine FrictionLaw (nvmax,keyopt,cinem,rvr,qa,deltat,qb,c,
     &                        intfer,kerr,hygro,ififo,founda,qab,ONED,
     &                        keyopt_1)
c --- primary function: calculate contact normal and tangent pressures by return mapping meth
c
c  input arguments:
c       nvmax   (int,sc,in)        - number of state variables
c       keyopt (int,ar(12),in)     - array of element key option
c       cinem  (dp,ar(14),in)      - kenimatic parameters used for contact
c         cinem(1) - structure moisture
c         cinem(2) - foundation moisture
c         cinem(3) - irigf = type of contact
c         cinem(5) - nocon = contact indicator
c                    nocon.ne.0--> no intersection between normal and foundation
c                         .eq.0-->    intersection between normal and foundation
c                         no contact at this integration point (nocon = 1 or 2)
c         cinem(6) - nocot = contact indicator at the step beginning
c         cinem(7) - indic = 0  mechanical,hygro stresses are computed
c                               compliance matix is computed
c         cinem(8) - kmst  = moisture key
c         cinem(9) - reserved
c         cinem(10)- frctn = frictional coefficient for isotropic friction
c                          = mean frictional coefficient for ortho. friction
c         cinem(11)-       - variational penatration factor
c         cinem(12)-       - penatration tolenent factor
c         cinem(13)-       - length of contact element
c         cinem(14)-       - element number
c         cinem(15)-       - pinb
c         cinem(16)-       - slave number
c         cinem(17)-       -       = segment element number
c                                  < slave element death
c         cinem(34)- frctnr= ratio of fric. coeff. 1 to mean fric. coeff.
c                             for ortho. friction
c         cinem(35)- frctns= ratio of fric. coeff. 2 to mean fric. coeff.
c                             for ortho. friction
c         cinem(36)- tnmax = maximum normal traction (normal bond strength)
c         cinem(37)- gnmax = maximum normal separation (gap at full normal
c                            damage/failure)
c         cinem(38)- ttmax = maximum tangential traction (tangetial bond strength)
c         cinem(39)- gtmax = maximum tangential separation (cum. slip at full
c                            normal damage/failure)
c         cinem(40)- vrpar = viscous regularization parameter for damage variable
c         cinem(41)- icflg = flag for tang. separation under normal compression
c                          = 0 - no tang. separation under compression (default)
c                          = 1 - tang. separation under normal compression with
c                                slip
c
c         rvr    (dp,ae(*),in) - element real constants
c
c         rvr(1)  - frigp  = penalty coefficient factor on contact pressure
c         rvr(2)  - ftoln  = tolence of penetration
c         rvr(3)  - igtol  = initial gap tolence
c         rvr(7)  - taumax = max value of contact friction
c         rvr(8)  - intf   = initial interference
c         rvr(9)  - opsf   = opening stiffness
c         rvr(10) - frigt  = tangent penalty coefficient factor
c         rvr(11) - b      = cohesion
c         rvr(16) - fact   = dynamic friction factor
c         rvr(24) - p      = adhesion
c         rvr(24) -        = scaling of interaction properties
c
c       qa     (dp,ar(*),in)       - stresses and state variables at initial
c                                    configuration
c       deltat (dp,sc,in)          - time increment
c
c  output arguments:
c       qb     (dp,ar(nmvax),out)  - stresses and state variables at updated
c                                    configuration
c        qb   (1) = presb  - normal contact pressure
c        qb   (2) = taurb  - tangent contact stress 1
c        qb   (3) = tausb  - tangent contact stress 2
c        qb   (5) = ield   -  1: if state is sticking contact (elastic)
c                          -  2: if state is sliding contact (elastoplastic)
c                          -  3: open but near contact
c                          -  4: open and not near contact
c        qb   (6) = ddisp  - energy dissipated per unit area due to friction
c        qb   (8) = iseg   - attached segment number calculated by CalPeni
c        qb   (9) = db     - total penetration distance calculated by CalPeni
c        qb  (10) = nocot  - calculated by CalPeni for initial situation
c        qb  (11) = dajct  - jocoba value used for hygro
c        qb  (12) = vtfrb  - 1st relative tangential velocity
c        qb  (13) = vtfsb  - 2nd relative tangential velocity
c        qb  (14) = gapi   - initial gap
c        qb  (15) = iconv  - convergent key
c        qb  (16) =        - total energy density
c        qb  (17) =        - sliding energy density
c        qb  (18) = rigp   - normal penetration coefficient
c        qb  (19) = toln   - penetration tolerance
c        qb  (20) = dba    - interpenetration variation
c        qb  (21) =        - augment Lagrange multilier for normal stress
c        qb  (22) =        - augment Lagrange multilier for 1st tangent stress
c        qb  (23) =        - augment Lagrange multilier for 2nd tangent stress
c        qb  (24) = rigt   - tangent penetration coefficient
c        qb  (25) =        - adapted normal penalty coefficient
c        qb  (26) =        - adapted tangent penalty coefficient
c        qb  (27) =        - normal presure at which slip started
c        qb  (28) =        - 0: contact never happen
c                          - 1: contact has happened
c        qb  (29) =        - accumulated sliding in 1st tangent direction
c        qb  (30) =        - accumulated sliding in 2nd tangent direction
c        qb  (31) =        - absoult accumulated sliding in 1st tangent direction
c        qb  (32) =        - absoult accumulated sliding in 2nd tangent direction
c        qb  (34) =        - moisture
c        qb(nvmax-13)=     - energy released during debonding in mode II
c        qb(nvmax-12)=     - energy released during debonding in mode I
c        qb(nvmax-11)=     - visc. regularized damage parameter for debonding
c        qb(nvmax-10)=     - failure mode for debonding
c                          - 0: no debonding
c                          - 1: mode I (normal mode)
c                          - 2: mode II (tangential mode)
c                          - 3: mixed mode (mode I + mode II)
c        qb(nvmax-9) =     - damage parameter for debonding
c        qb(nvmax-8) =     - debonding fraction, time
c        qb(nvmax-7) =     - initial interference
c        qb(nvmax-6) =     - old status
c        qb(nvmax-5) =     - old gap
c        qb(nvmax-4) =     - old normal pressure
c        qb(nvmax-3) =     - old 1st tangent velocity
c        qb(nvmax-2) =     - old 2nd tangent velocity
c        qb(nvmax-1) =     - old 1st friction stress
c        qb(nvmax) =       - old 2nd friction stress
c
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
#include "conv_stat.inc"
#include "dyn_info.inc"
c
      integer nvmax,keyopt(12),ififo(*)
c
      double precision cinem(43),rvr(W_NumContParam),
     &                 qa(nvmax),qb(nvmax),deltat, c(4,4),
     &                 dperr(10),hygro(16),founda(*), qab(*)

      integer mlaugm,iconv,nocon,ield,ield0,
     &        icut,aupp,irigf,kerr,irough,itie,knup,intfer,
     &        keyopt_1,ield1,icflg,idbmode,isflg
c
      double precision ZERO,ONE,HALF,TWO,rigp,rigt,b,taumax,p,
     &  toln,dba,frctn,da,db,db0,mdb,
     &  presa,taura,taurb,ddisp,presb0,ftol,
     &  presb,fplas,frigp,frigt,PREC,
     &  tausa,tausb,vtfra,vtfsa,vtfrb,vtfsb,
     &  rigtr,dar,dsr,fpres,dlam,tauab,taufact,rlh,fact,
     &  rigp0,rigt0,opsf,delueq,sdtol,fks,rigtnew,frctn0,closedgap,
     &  cohef,delasl,drigp,tnmax,ttmax,presc0,taurc0,tausc0,dfunc,
     &  dener1,dener2,dvnew,tauaa,sfac,gnmax,gtmax,crsher,crnor,
     &  ONED
c
      parameter (ZERO=0.0d0,ONE=1.0d0,HALF=0.5d0,TWO=2.0d0,PREC=1.0d-1)

#if defined(DEBUG)
      call SubBeg('FrictionLaw')
#endif

c --- initialize
c
      presc0 = ZERO
      taurc0 = ZERO
      tausc0 = ZERO

      idbmode = 0
      dfunc = ZERO
      dener1 = ZERO
      dener2 = ZERO
      dvnew = ZERO
      tauaa = ZERO
      sfac = ZERO
      isflg = 0

c --- material parameters
c
      frigp  = rvr(1)
      if (abs(frigp).le.TINY) frigp = ONE
      frigt = frigp
      toln = rvr(2)
      taumax = rvr(7)
      opsf = rvr(9)
c --- default estimation of opening stiffness
c      if (opsf.eq.ZERO) opsf = -qab(5)*qab(2)*0.2d0
      b = abs(rvr(11))
      cohef  = rvr(13)
      p = abs(rvr(24))
      fact   = abs(rvr(19))
c --- scale both adhesion and cohesion. same is done with debonding
      if (rvr(25).ne.ZERO) then
          p = abs(rvr(25))*p
          b = abs(rvr(25))*b
      endif
      if (fact.lt.ONE) fact = ONE
      if (A_KeyFirstPass.eq.1) then
          cinem(22) = ZERO
      endif

c --- cohesive zone model data for debonding
c
      tnmax = cinem(36)
c --- if tnmax = 0, reset tnmax to a small negative value
c     for mode II failure
      if (abs(tnmax).le.TINY) then
          tnmax = -sqrt(TINY)
      endif
      gnmax = cinem(37)
      ttmax = cinem(38)
      gtmax = cinem(39)
      icflg = nint(cinem(41))
      if (gnmax.le.ZERO .and. gtmax.le.ZERO) then
          icflg = 2
      elseif (gnmax.gt.ZERO .and. gtmax.gt.ZERO .and.
     &                            qa(nvmax-9).lt.ONED) then
c --- compute modification factor for separation distance
c     for debonding with mixed mode failure
          gnmax = gtmax/(ttmax/tnmax)
          cinem(37) = gnmax
          isflg = 0
      endif

c --- geometry parameters
c
      irigf = nint(cinem(3))
      nocon = nint(cinem(5))
      frctn = cinem(10)

c --- check for decohesion condition
c     if decohesion is met, friction is activated
c     otherwise not friction positive value for stress
c     negative value for total sliding distance
c     this is used with cohesion only
      if (cohef .lt. -TINY) then
          delasl = sqrt((qa(29) + qb(12))*(qa(29) + qb(12))
     &                 +(qa(30) + qb(13))*(qa(30) + qb(13)))
          if(delasl .lt. abs(cohef)) then
             frctn = ZERO
          endif
      endif
c
      frctn0 = frctn

c --- rigid friction law
c
      if (abs(frctn).le.TINY .and. b.ge.TINY) then
          frctn = 1.d-9
          frctn0 = 1.d0
      endif

c      rlh = cinem(12)/cinem(13)
      dperr(1) = cinem(14)
      call v_zero (hygro(1),16)
c
      aupp = keyopt(3)
      
      if (aupp.ge.3) aupp = aupp - 3
      if (A_KeyMatUpd.eq.1) aupp = 0
      db = qb(9)
      da = qa(9)
      db0 = db
      vtfra = qa(12)
      vtfrb = qb(12)
      vtfsa = qa(13)
      vtfsb = qb(13)
      if (nint(qa(5)).eq.0) qa(5) = 4.0d0
      ield0 = nint(qa(5))
      if (A_IterNum.eq.1) then
          ield1 = abs(ield0)
      else
          ield1 = nint(abs(qb(5)))
          if (ield1.eq.0) ield1 = 4
      endif
      if (ield0.eq.0) ield0 = 4
c
      presa = qa(1)
      taura = qa(2)
      tausa = qa(3)
c
      irough = 0
      itie = 0
      keyopt_1 = keyopt(1)

c --- if bonded surface has debonded completely,
c     set keyopt(1) locally to 0 (standard contact)
c
      if (abs(qa(nvmax-9)).ge.ONED .and. keyopt_1.ne.2) keyopt_1 = 0

c --- rough or bonded contact
c
      if (keyopt_1.eq.1.or.keyopt_1.eq.3.or.keyopt_1.ge.5) then
          irough = 1
          frctn = 1.d0
          frctn0 = 1.d0
          cinem(10) = 1.d0
      endif

c --- closed contact in first iteration if "no separation" or "bonded" contact or
c     "bonded" contact and small gap
c
      if (A_KeyFirstPass.eq.1) then
          qa(28) = ZERO
          if (keyopt_1.eq.4.or.keyopt_1.eq.5) qa(28) = ONE
          if (keyopt_1.eq.6.and.db.lt.1.d-4*cinem(13)) qa(28) = ONE
      endif
      qb(28) = qa(28)

c --- no separation or bonded contact
c
      if (qa(28).gt.HALF .and. keyopt_1.ge.2) itie = 1

c --- set/compute maximum penetration tolerance
c
      if (abs(toln).le.TINY) then
          toln = cinem(12)*PREC
          toln = min(toln,cinem(12),cinem(13),cinem(15))
      else
          if (toln.gt.ZERO) toln = cinem(12)*toln
          toln = min(toln,cinem(12),cinem(13),cinem(15))
          if (toln.lt.ZERO) toln = -toln
      endif

c --- set taumax
c
      if (taumax.le.TINY) taumax = 1.d20
      if (cinem(10).le.TINY .and. abs(taumax).ne.b)
     &    taumax = b

c --- set allowable number of augmentations
c
      mlaugm = 0
      if (keyopt(2).eq.1) then
          mlaugm = 4
      endif

c --- output C_nAugm (number of augmentations)
c
      if (A_IterNum.eq.1) then
          qa(27) = ZERO
          qb(27) = ZERO
      else
          if (C_nAugm.ne.0) qb(27) = C_nAugm
      endif

c --- initial configuration--open contact, update the penetration coefficient
c
      knup = 0

c --- set flag to update normal penalty stiffness
c --- first iteration of first substep
c
      if (A_IterNum.eq.1.and.A_SubStepNum.eq.1) knup = 1

c --- contact stiffness update
c
      if (aupp.ge.1) knup = 1

c --- open contact or ZERO normal pressure or small normal penalty stiffness
c
      if (nint(qa(5)).ge.3.or.qa(1).eq.ZERO.or.qa(18).le.TINY) knup = 1

c --- if history variables are updated,
c --- do not update normal penalty stiffness
c
      if (A_KeyMatUpd.eq.1) knup = 0

c --- compute normal penalty stiffness
c
      if (frigp.ge.0) then
          rigp = frigp*cinem(11)
          if (aupp.ge.1) then
              if (keyopt(6).eq.1) then
                  rigp = rigp*sqrt(cinem(12)*PREC/toln)
              elseif (keyopt(6).ge.2) then
                  rigp = rigp*cinem(12)*PREC/toln
              endif
          endif
      else
          rigp = - frigp
      endif
c
      if (knup.eq.1) then

c --- update normal penalty stiffness
c
c --- compute tangential penalty stiffness
c
         rigt = frctn0*cinem(11)
         taufact = frigt !*rlh
         if (rvr(10).gt.ZERO) taufact = rvr(10)*taufact
         if (taumax.gt.1.d-8) rigt = min(rigt,taumax/toln)
         if (taufact.le.ZERO) then
             rigt = -taufact*frctn0
         else
             rigt = taufact*rigt
         endif
         rigt = max(0.005d0*rigp,rigt)
         rigt = min(10.0d0*rigp,rigt)
         if (rvr(10).lt.ZERO) rigt = -rvr(10)

c         write (8,*) 'rigp, rigt',rigp, rigt, ififo(12), cinem(11)
         
c --- updated penalty method
c
         if (aupp.ge.3) then
             if (A_KeyFirstPass.eq.1.or.qa(18).le.TINY) then
                 qa(18) = rigp
                 qa(24) = rigt
             else
                 if (A_SubStepNum.eq.1.and.A_IterNum.eq.1) then
                     if (keyopt(3).ge.5) then
                         qa(18) = rigp*qa(25)
                         qa(24) = rigt*qa(27)
                     else
                         qa(18) = rigp*founda(25)
                         qa(24) = rigt*founda(27)
                     endif
                 else
                     if (keyopt(3).ge.5) then
                         qa(18) = qa(18)*qa(25)
                         qa(24) = qa(24)*qa(27)
                     else
                         qa(18) = qa(18)*founda(25)
                         qa(24) = qa(24)*founda(27)
                     endif
                     qa(18) = min(qa(18),1000.0d0*rigp)
                     qa(18) = max(qa(18),0.001d0*rigp)
                     qa(24) = min(qa(24),1000.0d0*rigt)
                     qa(24) = max(qa(24),0.001d0*rigt)
                 endif
c                 write (8,*) 'Update'
                 rigp = qa(18)
                 rigt = qa(24)
             endif
         endif
      else
c --- do not update normal penalty stiffness
c
         if (abs(qb(5)).le.2) then
             if (qa(18).gt.TINY) rigp = qa(18)
         endif
         rigt = ZERO
         if (qa(24).gt.TINY) rigt = qa(24)
      endif

      if (rvr(1).lt.-TINY)  then
          rigp = -rvr(1)
          if (tnmax.gt.TINY .and. db.gt.ZERO) then
              if (cinem(42).gt.ZERO) then
                  mdb = HALF*rigp/cinem(42)
                  mdb = min(mdb,db)
                  rigp = rigp - mdb*cinem(42)
              elseif (cinem(42).lt.ZERO) then
                  rigp = rigp + db*cinem(42)
              endif
          endif
      endif
      if (rvr(10).lt.-TINY) then
          rigt = -rvr(10)
          if (ttmax.gt.TINY) then
              delueq = sqrt(qb(29)*qb(29)+qb(30)*qb(30))
              if (cinem(43).gt.ZERO) then
                  mdb = HALF*rigt/cinem(43)
                  mdb = min(mdb,delueq)
                  rigt = rigt - mdb*cinem(43)
              elseif (cinem(43).lt.ZERO) then
                  rigt = rigt + delueq*cinem(43)
              endif
          endif
      endif
      
c --- generalized stresses time integration
c
c
c --- compute total and increm. penetration and tangential sliding
c
      if (db.ge.ZERO) then
          if (db.lt.sqrt(TINY)*cinem(13)) db = -sqrt(TINY)*cinem(13)
          if (A_KeyFirstPass.eq.1.and.keyopt(8).eq.0) then
             if (qb(nvmax-7).lt.ZERO.and.db.lt.1.d-7*cinem(13))
     &           db = -sqrt(TINY)*cinem(13)
         endif
      endif

      dba  = -(db-da)
c --- calculate mechanical ! and hygro stresses - future work
c
      ddisp = ZERO
      presb = ZERO
      taurb = ZERO
      tausb = ZERO
      tauab = ZERO
      rigtr = ZERO
      dar = ZERO
      dsr = ZERO
      delueq = sqrt(vtfrb*vtfrb+vtfsb*vtfsb)
      sdtol = 0.01d0*cinem(13)

c --- store normal and tangential penalty stiffnesses
c
      qa(18) = rigp
      qa(24) = rigt
c
      rigp0 = rigp
      if (isflg.eq.0) then
          rigt0 = rigt
      else
          rigt0 = sfac*rigp
      endif
      
c --- compute friction coefficient based on exponential decay
c
      if (rvr(19).gt.ONE) then
          if (abs(deltat).gt.sqrt(TINY)) then
              frctn =frctn*(ONE+(fact-ONE)*exp(-rvr(20)*delueq/deltat))
          endif
          if (rvr(20) .eq. ZERO .and. ield0 .eq. 2) then
              frctn = frctn / fact
          endif
      endif
      cinem(10) = frctn
      ield  = 4

c --- set initial status
      if (A_KeyFirstPass.eq.1)  then
         if (db.gt.ZERO.or.nocon.ne.0) then
c --- gap open
            if (A_Debug.ge.2) then
               dperr(2) = db
               call erhandler('FrictionLaw',5000,1,'Contact status of
     &              element %I is initially open. The gap is %g.',
     &              dperr(1),' ')
            endif
            if (irigf.eq.0.and.abs(qb(9)).lt.cinem(15)) ield = 3
            if (nocon.eq.0) then
                ield = 3
                if (abs(qb(9)).lt.cinem(15)) then
                    if (keyopt_1.eq.4) ield = 2
                    if (keyopt_1.eq.5) ield = 1
	          endif
            endif
         else
c --- gap closed
            icut = 0
            if (irigf.ge.1.and.-db.ge.HALF*toln) icut = 1
            if (irigf.ge.2.and.-db.ge.0.1d0*toln) icut = 1
            if (icut.eq.0) then
                ield = 2
                if (keyopt_1.eq.5) ield = 1
                if (A_Debug.ge.2) then
                    dperr(2) = db
                    call erhandler('FrictionLaw',5000,1,'Contact status
     &                             of element %I is initially closed.
     &                             The penetration is %g.',dperr(1),' ')
                endif
                if (delueq.le.sqrt(TINY).and.frctn.ge.TINY) ield = 1
            endif
         endif
         qa(5) = ield
         qb(5) = ield
         ield0 = ield
      endif
      ield  = 4

      if (irigf.eq.0.and.abs(qb(9)).lt.cinem(15)) ield = 3
     
      if (nocon.eq.0) then
          ield = 3
          if (abs(qb(9)).lt.cinem(15)) then
       	      if (keyopt(1).eq.4) ield = 2
       	      if (keyopt(1).eq.5) ield = 1
          endif
      endif
c
      if ((keyopt(2).eq.0).and.(presa.eq.ZERO.and.db.ge.ZERO).and.
     &     itie.eq.0) goto 300
      if (A_KeyFirstPass.ne.1.and.keyopt_1.eq.6.and.
     &    qa(28).lt.HALF) goto 300

c --- augument stresses
c
      if (mlaugm.ge.1) then
         if (A_IterNum.eq.1) then
            qb(22) = qa(2)
            qb(23) = qa(3)
            qa(15) = 0
            qb(15) = 0
         else
            if (nint(qb(15)).eq.1.and.C_Disp*C_Forc.eq.1) then
               qa(15) = qa(15) + 1
            endif
         endif
         presa = qa(1)
         taura = qb(22)
         tausa = qb(23)
      endif
      if (nocon.ne.0) goto 300
      if (itie.eq.1.and.db.gt.ZERO) then
c --- check for debonding; do not change penalty stiffness
c --- if debonding is active with bonded contact 
cc --- unless it is a first pass
c
c         if (.not.(abs(tnmax).le.TINY) .and.
c     &       (keyopt_1.ge.4.and.keyopt_1.le.6)) then
c --- stiffness if openning contact
            if (opsf.lt.ZERO) then
                rigp0 = - opsf
            elseif (opsf.gt.TINY) then
                rigp0 = opsf*rigp
            endif
c         endif
      endif

c --- elastic contact : no sliding
c
      if (mlaugm.ge.1) then
          if (-db.gt.toln.and.A_IterNum.gt.1.and.A_KeyMatUpd.ne.1)then
               qb(21) = qb(21) - HALF*rigp0*db
          endif
          if (db.le.ZERO) then
              presb = qb(21) - rigp0*db
          else
              qb(21) = ZERO
              presb = - rigp0*db
          endif
      else
          qb(21) = ZERO
          presb = -rigp0*db
      endif

c --- maximum tensile force
      ftol = HUGE
      if (rvr(22).gt.TINY) then
          ftol = rigp0*rvr(22)
      elseif (rvr(22).lt.-TINY) then
          ftol = -rvr(22)
      endif

      if (-presb.gt.ftol) then
          presb = -ftol
      endif

      qb(nvmax-4) = rigp0

c --- calculate tangent stiffness
c
      rigtnew = rigt
      presb0 = presa

      if (aupp.eq.2) presb0 = presb
      fks = ONE
      if (-rvr(10).gt.TINY) then
          rigtnew = -rvr(10)
          goto 144
      endif

      if (keyopt_1.eq.1.or.keyopt_1.eq.3.or.keyopt_1.ge.5)
     &    sdtol = HALF*sdtol
      if (rvr(21).gt.TINY) then
          sdtol = rvr(21)*cinem(13)
      elseif (rvr(21).lt.-TINY) then
          sdtol = -rvr(21)
      endif
      if (aupp.eq.0) then
          rigtnew = rigt
          goto 144
      endif
      fks = rvr(10)
      if (fks.le.TINY) fks = ONE
      if (frctn0*abs(presb0).lt.taumax) then
          rigtnew = fks*frctn0*abs(presb0)/sdtol
      else
          rigtnew = fks*taumax/sdtol
      endif

      if (rigtnew.le.rigt*0.1d0) then
          rigtnew = 0.1d0*rigt
      endif
 144  continue

      if (knup.eq.1) then
          rigt = rigtnew
          qa(24) = rigt
      endif


      if (isflg.eq.0) then
          rigt0 = rigt
      else
          rigt0 = sfac*rigp
      endif

c --- DO NOT CHANGE TANGENT
c      if (itie.eq.1.and.db.gt.ZERO) then
c --- check for debonding; do not change penalty stiffness
c --- if debonding is active with bonded contact 
cc --- unless it is a first pass
c
c         if (.not.(abs(tnmax).le.TINY) .and.
c     &       (keyopt_1.ge.4.and.keyopt_1.le.6)) then
c --- stiffness if openning contact
c            if (opsf.lt.ZERO) then
c               rigt0 = - opsf
c            elseif (opsf.gt.TINY) then
c               rigt0 = opsf*rigt
c            endif
c         endif
c      endif

      if (itie.eq.1 .and. keyopt_1.ne.2) goto 150
c --- if bonded surface has debonded completely, change
c     contact status from tied to sliding contact,
c     or open contact (if debonding has occured)
c
      if (qa(nvmax-9).ge.ONED) then
          if (presb.gt.ZERO .and. keyopt_1.ne.2) then ! compression--contact
              if (icflg.eq.1) then
                  ield = 2          ! sliding
              endif
          else                    ! tension--separation in normal direction
              if (abs(qb(9)).lt.cinem(15)) then
                  ield = 3          ! open but near
              else
                  ield = 4          ! open and far
              endif
          endif

          if (keyopt_1.eq.2) then
c --- zero-out debonding history, start new
              presb = ZERO
              taurb = ZERO
              tausb = ZERO
              frctn = ZERO

              qb(21) = ZERO
              qb(28) = ZERO
              qb(nvmax-8) = ZERO
              qb(nvmax-9) = ZERO
              qb(nvmax-10) = ZERO
              qb(nvmax-11) = ZERO

              ield = 3          ! open but near
              goto 300
          endif
      endif

      if (itie.eq.1) goto 150

c --- save presb
      presb0 = presb
      presb = presb + p

      if (presb.le.ZERO) then
          if (presb.eq.ZERO.and.b.gt.ZERO) goto 145
          presb = ZERO
          goto 300
      endif

 145  continue
      if (presb.gt.ZERO) then
          fplas = frctn*presb + b
      else
          if (rvr(11).ge.ZERO) then
              fplas = b
          else
              fplas = max(ZERO,frctn*presb+b)
          endif    
      endif
c --- restore presb
      presb = presb0
      if (fplas.le.sqrt(TINY)*rigp0*cinem(13)) then
          ield  = 2
          if (frctn.gt.1.d-8) ield = 1
          goto 300
      endif

 150  continue

c --- special case: free sliding
c
      if (frctn.eq.ZERO .and. taumax.eq.ZERO) then
          ield  = 2
          goto 300
      endif

      if ((keyopt_1.eq.2.or.keyopt_1.eq.4).and.(presb+p).lt.ZERO.and.
     &     b.le.ZERO) then
          ield  = 2
          frctn = ZERO
          goto 300
      endif
      ield  = 1
c
      presb0 = presb + p
      
      if (ield.eq.1) then
          taurb = (qa(29) + qb(12))*rigt0
          tausb = (qa(30) + qb(13))*rigt0
      else
          taurb = taura + vtfrb*rigt0
          tausb = tausa + vtfsb*rigt0
      endif

      qb(nvmax-3) = rigt0
c
      tauab = sqrt(taurb*taurb+tausb*tausb)
      if (irough.eq.1) then
          if (abs(deltat).gt.sqrt(TINY))
     &    ddisp = abs(tauab*delueq)/deltat
          goto 300
      endif
      if (tauab.lt.abs(presb)*1.d-8.and.frctn.lt.TINY) then
          ield = 2
          goto 300
      endif
c
      if (presb0.gt.ZERO) then
          fpres = frctn*presb0 + b
      else
          if (rvr(11).ge.ZERO) then
              fpres = b
          else
              fpres = max(ZERO,frctn*presb0+b)
          endif    

      endif

c --- check if fpres => maximum imposed contact friction
c
      if (fpres.ge.taumax ) then
          fpres = taumax
      endif
      
c --- check whether plastic deformation occurs
c
      fplas = tauab - fpres
c
      if (fplas.ge.0) goto 200

c --- elastic state
      ield = 1
      goto 300

c --- plastic contact : sliding
c
 200  continue
c
      ield = 2
      if (abs(tauab).le.TINY.and.frctn.gt.TINY) then
          if (frctn.gt.TINY) then
              ield = 1
          else
              ddisp = ZERO
              taurb = ZERO
              tausb = ZERO
              tauab = ZERO
              rigtr = ZERO
              dar = ZERO
              dsr = ZERO
          endif
          goto 300
      endif
      dlam  = fpres / tauab
c
      dar   = taurb / tauab
      dsr   = tausb / tauab
c
      taurb = taurb * dlam
      tausb = tausb * dlam
      tauab = sqrt (taurb*taurb + tausb*tausb)
c
      rigtr = rigt0  * dlam
c
      if (abs(deltat).gt.sqrt(TINY)) ddisp = abs(tauab*delueq)/deltat
c
  300 continue
c
      if (ield.eq.4) then
          if (nocon.eq.0) ield = 3
          if (irigf.eq.0.and.abs(qb(9)).lt.cinem(15)) ield = 3
      endif
      if (presb.lt.ZERO) then
          if (ield.eq.1) ield = -1
          if (ield.eq.2) ield = -2
      endif
      qb(5) = abs(ield)
      if (abs(ield).gt.2) qb(21) = ZERO

c --- case of self spurious contact
      if (irigf.ne.0.and.ield.le.2.and.ield0.eq.4.and.
     &    -db.gt.0.2d0*cinem(15)) then
c          write (6,*) 'Leg S0', ield, ield0, ififo(12), cinem(15)
          if (C_ConnUpd.eq.0.or.kerr.eq.3) then
            if (C_nConF.eq.0) then
               call lockSMP (W_LockConvVar)
               C_nConF = C_nConF + 1
               call unlockSMP (W_LockConvVar)
               if (A_Debug.eq.1) then
                 dperr(1) = cinem(14)
                 dperr(2) = cinem(16)
                 call erhandler ('FrictionLaw',5000,2,'Contact element
     &               %I status changes suddenly with slave element %I.
     &               It is a spurious contact.',dperr(1),' ')
               endif
            endif
          endif
          ddisp = ZERO
          presb = ZERO
          qb(21) = ZERO
          taurb = ZERO
          tausb = ZERO
          rigtr = ZERO
          dar = ZERO
          dsr = ZERO
          ield = 4
          if (A_SpecVer.ne.1) then
              qb(5) = 5.0d0
          else
              qb(5) = 4.0d0
          endif
          qb(8) = 0
          qb(9) = cinem(15)
          db0 = cinem(15)
      endif
      if (irigf.ne.0.and.ield.le.2.and.ield0.eq.3.and.
     &    -db.gt.HALF*cinem(15)) then
c          write (6,*) 'Leg S1', ield, ield0, ififo(12), cinem(15), db
          if (C_ConnUpd.eq.0.or.kerr.eq.3) then
             if (C_nPen.eq.0) then
                 call lockSMP (W_LockConvVar)
                 C_nPen = C_nPen + 1
                 call unlockSMP (W_LockConvVar)
                 if (A_Debug.eq.1) then
                     dperr(1) = cinem(14)
                     dperr(2) = -db
                     dperr(3) = cinem(16)
                     call erhandler ('FrictionLaw',5000,2,
     &              'Contact element %i has too much penetration %g
     &               related to slave element %i. It may be a
     &               spurious contact.',dperr(1),' ')
                 endif
             endif
          endif
          ddisp = ZERO
          presb = ZERO
          qb(21) = ZERO
          taurb = ZERO
          tausb = ZERO
          rigtr = ZERO
          dar = ZERO
          dsr = ZERO
          ield = 4
          if (A_SpecVer.ne.1) then
              qb(5) = 5.0d0
          else
              qb(5) = 4.0d0
          endif
          qb(8) = 0
          qb(9) = cinem(15)
          db0 = cinem(15)
      endif
      if (intfer.eq.1) goto 350
      if ((ield.eq.1.or.ield.eq.2) .and. irigf.ne.0 .and.
     &    A_KeyFirstPass.eq.1) then
         icut = 0
         if (irigf.ge.1.and.-db0.ge.0.1d0*cinem(15)) icut = 1
         if (irigf.ge.2.and.-db0.ge.HALF*toln) icut = 1
         if (icut.eq.1) then
             ddisp = ZERO
             presb = ZERO
             qb(21) = ZERO
             taurb = ZERO
             tausb = ZERO
             rigtr = ZERO
             dar = ZERO
             dsr = ZERO
             ield = 4
             qb(8) = 0
             qb(9) = cinem(15)
             if (A_SpecVer.ne.1) then
                 qb(5) = 5.0d0
             else
                 qb(5) = 4.0d0
             endif
             if (C_nPen.eq.0.and.(C_ConnUpd.eq.0.or.kerr.eq.3)) then
                 call lockSMP (W_LockConvVar)
                 C_nPen = C_nPen + 1
                 call unlockSMP (W_LockConvVar)
                 if (A_Debug.eq.1)  then
                     dperr(1) = cinem(14)
                     dperr(2) = -db0
                     dperr(3) = cinem(16)
                     call erhandler ('FrictionLaw',5000,2,
     &                 'Contact element %i has too much penetration %g
     &                  related to slave element %i. It may be a
     &                  spurious contact.',dperr(1),' ')
                 endif
            endif
            db0 = cinem(15)
         endif
      endif
 350  continue

      qb(1) = presb
      qb(2) = taurb
      qb(3) = tausb
      tauab = sqrt (taurb*taurb + tausb*tausb)
      qb(4) = 0
      qb(5) = abs(ield)
      if (C_ConnUpd.eq.0) then
          if ((abs(ield).le.2 .and. abs(ield1).ge.3) .or. 
     &        (abs(ield1).le.2 .and. abs(ield).ge.3)) then
c             qb(44) = qb(44) + abs(abs(ield)-abs(ield1))
c              write (8,*) 'Leg 3', ififo(12), qab(19),qab(20),qab(21)
          endif   
      endif
      qb(5) = abs(ield)
      qb(6) = ddisp
      qb(18) = qa(18)
      qb(19) = toln
      qb(20) = dba
      qb(24) = qa(24)
      if (abs(ield).ne.1.and.abs(ield).ne.2) then
          delueq = ZERO
          qb(29) = qa(29)
          qb(30) = qa(30)
          qb(31) = qa(31)
          qb(32) = qa(32)
      else
          qb(29) = qa(29) + qb(12)
          qb(30) = qa(30) + qb(13)
          qb(31) = qa(31) + abs(qb(12))
          qb(32) = qa(32) + abs(qb(13))
      endif
      if (abs(ield).le.2.or.(abs(ield).eq.3.and.qb(9).ge.-TINY)) then
c --- get the closest contact point information
         if (qb(9).lt.founda(9).and.cinem(16).gt.0) then
             call lockSMP (W_LockConvVar)
             founda(9) = qb(9)
             ififo(16) = cinem(14)
             ififo(17) = cinem(16)
             call unlockSMP (W_LockConvVar)
         endif
         closedgap = qb(9)+qb(nvmax-7)+rvr(8)
         if (keyopt(5).gt.0.and.keyopt(5).lt.4)
     &       closedgap = closedgap + founda(22)
         if (abs(ield).le.2.and.cinem(16).gt.0) then
            if (closedgap.gt.founda(24).and.closedgap.gt.ZERO) then
               call lockSMP (W_LockConvVar)
               founda(24) = closedgap
               if (abs(founda(24)).ge.abs(founda(23))) then
                   founda(25) = cinem(14)
                   founda(26) = cinem(16)
               endif
               call unlockSMP (W_LockConvVar)
            endif
            if (closedgap.lt.founda(23).and.closedgap.lt.ZERO) then
               call lockSMP (W_LockConvVar)
               founda(23) = closedgap
               if (abs(founda(23)).gt.abs(founda(24))) then
                   founda(25) = cinem(14)
                   founda(26) = cinem(16)
               endif
               call unlockSMP (W_LockConvVar)
            endif
         endif
      endif

      if (nint(qa(27)).eq.0.and.ield.eq.2) qa(27) = ONE

c --- set convergence key if contact with augmented Lagrange method
c
      iconv = 0
      if (ield.ge.3 .or. mlaugm.eq.0) goto 400
      if (-db0.gt.toln) iconv = 1
      if (db0.gt.toln.and.itie.eq.0) iconv = 1
c
  400 continue
      qb(15) = iconv

  410 continue
      qb(25) = qa(25)
      qb(26) = qa(26)
      if (abs(ield).eq.1.or.abs(ield).eq.2) then
          qb(28) = 1.0d0
          if (keyopt_1.ge.4.or.A_KeyFirstPass.eq.1) qa(28) = 1.0d0
      endif
      icut = 0

c --- check debonding
c
      if (keyopt_1.eq.0) goto 420
      if (qb(5).ne.1.and.qb(5).ne.2) goto 420
      if (abs(tnmax).gt.TINY) then
         if (icflg.ne.2) then
             idbmode = nint(qa(nvmax-10))
             if (idbmode.eq.0) then
                 if (ttmax.gt.TINY) then
                     if (tnmax.lt.TINY) then
                         idbmode = 2   !tnmax < 0 and ttmax > 0
                     else
                         idbmode = 3   !tnmax > 0 and ttmax > 0
                     endif
                 endif
                 if (tnmax.gt.TINY) then ! idbmode = 0 or 3, tnmax > 0
                     if (idbmode.eq.0) then
                         idbmode = 1
                     endif
                 endif
                 qb(nvmax-10) = idbmode
                 qa(nvmax-10) = qb(nvmax-10)
             endif

             if (idbmode.gt.0) then
                 call ContDamage(nvmax,cinem(1),qa(1),qb(1),kerr,
     &                           presb,taurb,
     &                           tausb,presc0,taurc0,tausc0,
     &                           rigp0,rigt0,idbmode,dfunc,ONED)
c --- update energy released during debonding
c
                 dvnew = qb(nvmax-11)

                 if (dvnew.gt.ZERO .and. dvnew.lt.ONED) then
                     tauaa = sqrt(taura*taura + tausa*tausa)
                     if (idbmode.eq.1) then
                         if (presb.lt.ZERO) then
                             dener1 = abs(HALF*(presa+presb)*dba)
                         else
                             dener1 = ZERO
                         endif
                     elseif (idbmode.eq.2) then
                         dener2 = abs(HALF*(tauaa+tauab)*delueq)
                     elseif (idbmode.eq.3) then
                         if (presb.lt.ZERO) then
                             dener1 = abs(HALF*(presa+presb)*dba)
                         else
                             dener1 = ZERO
                         endif
                         dener2 = abs(HALF*(tauaa+tauab)*delueq)
                     endif
                     qb(nvmax-12) = qa(nvmax-12) + dener1
                     qb(nvmax-13) = qa(nvmax-13) + dener2
                 endif
             endif
         else
c --- simplified mode
             if (ttmax.lt.TINY .or. keyopt_1.eq.2) then
                 crsher = ZERO
             else
                 crsher = tauab/ttmax
             endif
             if (tnmax.lt.TINY .or. presb.gt.ZERO) then
                 crnor = ZERO
             else
                 crnor = presb/tnmax
             endif
             crnor = sqrt(crnor*crnor+crsher*crsher)
             if (crnor.gt.ONE .and. A_KeyMatUpd.eq.1) then
                 presb = ZERO
                 tauab = ZERO
                 qb(nvmax-9) = crnor
                 qb(nvmax-6) = A_TimeCur
             endif
         endif

      endif

c --- energy
c
c --- disipated
      qb(17) = qa(17) + ddisp*deltat
c --- elastic
      qb(16) = ZERO
      if (abs(ield).gt.2) goto 410
      if (abs(ield).eq.2.and.rigp0.gt.1.d-10) qb(16) = qb(16) +
     &    HALF*presb*presb/rigp0
      if (abs(ield).eq.1.and.rigt0.gt.1.d-10) qb(16) = qb(16) +
     &    HALF*(presb*presb/rigp0 + tauab*tauab/rigt0)

 420  continue

c --- check for too big penetration
c
      if (intfer.eq.1) then
         if (A_SubStepNum.eq.1) goto 500
         goto 450
      endif
      if (ield.eq.1.or.ield.eq.2.and.C_ConnUpd.eq.0) then
         if (-db0.gt.30.d0*toln) icut = 1
c         if (-db0.gt.3.0d0*cinem(12)) icut = 1
         if (mlaugm.ge.1) then
             if (-db0.gt.8.0d0*toln.and.qa(15).ge.1) icut = 1
         endif
         if (icut.eq.1) then
c             write (6,*) 'Leg 00', db, 30*toln, 3.0d0*cinem(12)  
c             if (A_SpecVer.gt.0 .and. A_KeyMatUpd.eq.0) then
c                 call PutErr(1)
c                 dperr(1) = ififo(12)
c                 call lockSMP (W_LockAnsys)
c                 call erhandler ('FrictionLaw',5000,2,
c     &            'Too much penetration at interface %I. 
c     &             Bisection triggered.',dperr(1),' ')
c                 call unlockSMP (W_LockAnsys)
c                 goto 999
c             endif
c             write (6,*) 'Too much penetration',db0,ififo(12)
             dperr(1) = cinem(14)
             dperr(2) = ififo(12)
             dperr(3) = -db0
             dperr(4) = cinem(16)
             ield = 3
             qb(1) = ZERO
             qb(2) = ZERO
             qb(3) = ZERO
             if (kerr.eq.3) then
                 call erhandler ('FrictionLaw',5000,kerr,
     &               'Contact element %i (ID %i) has too much
     &                penetration %g related to slave element %i.
     &                Substep will be abandoned.',dperr(1),' ')
             elseif(C_nPen.eq.0.and.(C_ConnUpd.eq.0.or.kerr.eq.3)) then
                 call lockSMP (W_LockConvVar)
                 C_nPen = C_nPen + 1
                 call unlockSMP (W_LockConvVar)
                 if (A_Debug.eq.1)  then
                     call erhandler ('FrictionLaw',5000,kerr,
     &                  'Contact element %i (ID %i) has too much
     &                   penetration %g related to slave element %i.',
     &                   dperr(1),' ')
                 endif
             endif
             goto 500
         endif
      endif
 450  continue

c --- check for changing contact status
c
      if (ield0.le.2.and.abs(ield).gt.2) then
c --- pair got disconnected
          if (C_ConnUpd.eq.0) then
              if (C_nConF.eq.0) then
                  call lockSMP (W_LockConvVar)
                  C_nConF = C_nConF + 1
                  call unlockSMP (W_LockConvVar)
              endif
          endif
      endif
      
      if (ield0.eq.4.and.abs(ield).le.2.and.
     &   (C_ConnUpd.eq.0.or.kerr.eq.3)) then

c         write (6,*) 'Leg special',cinem(14),cinem(16)
         ield = 4  ! mark for deletion
         qb(1) = ZERO
         qb(2) = ZERO
         qb(3) = ZERO
         if (A_SpecVer.ne.1) then
             qb(5) = 5.0d0
         else
             qb(5) = 4.0d0
         endif
             
         qb(44) = qb(44) + 2
c         if (A_SpecVer.gt.0 .and. A_KeyMatUpd.eq.0) then
c             call PutErr(1)
c             dperr(1) = ififo(12)
c             call lockSMP (W_LockAnsys)
c             call erhandler ('FrictionLaw',5000,2,
c     &           'Sudden contact status change at interface %I.
c     &            Bisection triggered.',dperr(1),' ')
c            call unlockSMP (W_LockAnsys)
c         endif
         if (C_nConF.eq.0) then
             call lockSMP (W_LockConvVar)
             C_nConF = C_nConF + 1
             call unlockSMP (W_LockConvVar)
             if (A_Debug.eq.1) then
                 dperr(1) = cinem(14)
                 dperr(2) = ififo(12)
                 dperr(3) = cinem(16)
                 call erhandler ('FrictionLaw',5000,kerr,'Master
     &           element %I (ID %I) changed status suddenly from
     &           no-contact to contact (with) slave element %I).',
     &           dperr(1),' ')
             endif
         endif
         goto 999
      elseif (ield0.le.2.and.ield.eq.4.and.cinem(17).gt.0.and.
     &       (C_ConnUpd.eq.0.or.kerr.eq.3)) then
         ield = 4  ! mark for deletion
         qb(1) = ZERO
         qb(2) = ZERO
         qb(3) = ZERO
         if (A_SpecVer.ne.1) then
             qb(5) = 5.0d0
         else
             qb(5) = 4.0d0
         endif
         qb(44) = qb(44) + 1
         if (C_nConF.eq.0) then
             call lockSMP (W_LockConvVar)
             C_nConF = C_nConF + 1
             call unlockSMP (W_LockConvVar)
             if (A_Debug.eq.1) then
                 dperr(1) = cinem(14)
                 dperr(2) = ififo(12)
                 dperr(3) = cinem(17)
                 call erhandler ('FrictionLaw',5000,kerr,'Master element
     &           %I (ID %I) changed status suddenly from contact (with
     &           slave element %I) to no-contact.',dperr(1),' ')
             endif
         endif
      endif
c
 500  continue
      qa(25) = ONE
      qa(26) = ONE

c --- select penalty coefficient for next iteration/step
      if (ield.le.2.and.aupp.ge.3) then
          if (A_KeyMatUpd.eq.1) then
c --- convergence has reached
c --- select penalty coefficient for next step
c --- increase penalty coefficient if too much penetration
             if (A_IterNum.lt.3) then
                 qa(25)=1.5d0
                 if (ield.eq.1) qa(26)=1.5d0
                 founda(24) = 1.5d0
                 if (ield.eq.1) founda(26) = 1.5d0
             elseif (A_IterNum.gt.4) then
c --- decrease penalty coefficient if bad converence
                 if (-qb(9).lt.toln) qa(25)=0.75d0
                 if (ield.eq.1.and.delueq.lt.sdtol) qa(26)=0.75d0
             endif
         else
c --- select penalty coefficient for next iteration
c --- increase penalty coefficient if too much penetration
             if (-qb(9).gt.toln) then
                 if (-qb(9).gt.-0.25d0*qa(nvmax-5)) then
                     qa(25) = TWO
                     founda(24) = TWO
                 endif
             endif
             if (ield.eq.1.and.delueq.gt.sdtol) then
                 if (delueq.gt.sqrt(qa(nvmax-2)*qa(nvmax-2)+
     &               qa(nvmax-3)*qa(nvmax-3))) then
                     qa(26) = TWO
                     founda(26) = TWO
                 endif
             endif
             qa(nvmax-5) = qb(9)
             qa(nvmax-2) = qb(12)
             qa(nvmax-3) = qb(13)
         endif
         if (founda(24).gt.qa(25)) then
             call lockSMP (W_LockConvVar)
             founda(24) = qa(25)
             call unlockSMP (W_LockConvVar)
         endif
         if (ield.eq.1.and.founda(26).gt.qa(26)) then
             call lockSMP (W_LockConvVar)
             founda(26) = qa(26)
             call unlockSMP (W_LockConvVar)
         endif
      endif
      qb(25) = qa(25)
      qb(26) = qa(26)

c --- calculate consistent material tensor
c
      call v_zero (c(1,1),16)

      if (abs(ield).eq.1.or.abs(ield).eq.2) goto 900
      if (keyopt_1.gt.1.or.rvr(12).ge.-TINY) goto 900
c --- add contact damping
c     assumes normal and tangential direction have the same contact stiffness
c      if (nocon.eq.0 .and. qa(nvmax-9).ne.-ONE) then
      if (nocon.eq.0) then
          if (db.ge.cinem(15)) goto 900
          if (abs(deltat).lt.sqrt(TINY)) goto 900 
          ield = -2
          qb(5) = dble(ield)
          rigp0 = -rvr(12)
          drigp = ZERO
          if (db.gt.HALF*cinem(15)) then
              rigp0 = TWO*rvr(12)*(db-cinem(15))
              drigp = TWO*rvr(12)/cinem(15)*db
          endif

          if (A_Trans.gt.0) then
              rigp0 = rigp0*I_Prm(2)
          else
              rigp0 = rigp0/deltat
          endif
         
          presb = -rigp0*db
          taurb = rigp0*vtfrb
          tausb = rigp0*vtfsb
          qb(1) = presb
          qb(2) = taurb
          qb(3) = tausb
          qb(21) = ZERO
          c(1,1) = rigp0
          c(1,3) = drigp*(vtfrb-vtfra)/deltat
          c(2,2) = rigp0
          c(2,3) = drigp*(vtfsb-vtfsa)/deltat
          c(3,3) = rigp0-drigp*dba/deltat
          goto 999
      endif

 900  continue

      if (abs(ield).ne.1.and.abs(ield).ne.2) goto 999

c --- calculate tangent stiffness for debonding
c
       if (idbmode.gt.0) then
           call ContDebond(nvmax,cinem(1),qa(1),qb(1),
     &                     presc0,taurc0,tausc0,dfunc,rigp0,
     &                     rigt0,rigt0,idbmode,c(1,1),ONED)
           goto 999
       endif

c --- definition of c: sigma = c * epsilon
c
      if (abs(ield).eq.2) then
          c(1,1) = rigtr * (ONE - dar*dar)
          c(2,2) = rigtr * (ONE - dsr*dsr) 
      else
          c(1,1) = rigt0
          c(2,2) = rigt0
      endif 
      c(3,3) = rigp0
      
 999  continue
      
c --- if defined debonding, we only count the first contact
      if (keyopt_1.eq.3.and.((cinem(36)+cinem(38)).gt.TINY)) then
          if (A_KeyConv.eq.1.and.A_KeySvrUpd.eq.1) then
              if (abs(ield).ne.1) then
                  qa(nvmax-9) = -ONE
                  qb(nvmax-9) = -ONE
              endif
          endif
      endif

#if defined(DEBUG)
      call SubEnd('FrictionLaw')
#endif

      return
      end
