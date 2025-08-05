      subroutine MaterialDriver (ireal, IntId, keycut, 
     &                           matRecIndex,  matRec, matMWIndex, 
     &                           matMW, matWorkIndex, matWork,
     &                           strainInc, moist, nTens, elemId, matId,
     &                           kDomIntPt, kSectPt, kMatPoint,
     &                           mref, coords, basis,
     &                           iniFlag, tMatProp, khygro,
     &                           stress, consJac, strain,
     &                           siglv, siglvi, eEnerPt, nlkword
#if defined(MOIST)
     &                           ,eh)
#else
     &                          )
#endif            
c --- primary function: this is the master material interface; given the state of
c                       the material at time 't' and incremental deformation
c                       update the state of the material at 't+Dt' and pass
c                       back the consistent material Jacobian
c
c  input arguments:
c       keycut      (int,sc,in)              - bisec/cut flag
c       matRecIndex (int,ar(*),in)           - array indexing into matRec
c       matRec      (dp,ar(*),in)            - material point saved variables
c       matMWIndex  (int,ar(*),in)           - array indexing into matMW array
c       matMW       (dp,ar(*),in)            - mat work array used to modify matRec at convergence
c       matWorkIndex(int,ar(*),in)           - array indexing into matWork array
c       matWork     (dp,ar(*),in)            - mat work array (never saved)
c       strainInc   (dp,ar(nTens),in)        - strain increment
c       moist       (dp,sc,in)               - moisture at the material point
c       nTens       (int,sc,in)              - nDirect + nShear
c       elemId      (int,sc,in)              - element number
c       matId       (int,sc,in)              - material number
c       kMatPoint   (int,sc,in)              - "k"th domain integration point
c       kSectPt     (int,sc,in)              - "k"th Section point
c       mref        (dp,sc,in)               - reference moisture
c       coords      (dp,ar(3),in)            - coordinates of this pt
c       basis       (dp,ar(3,3),in)          - basis directions at this pt
c       iniFlag     (log,sc,in)              - initial stress flag
c
c  output arguments:
c       tMatProp    (dp,ar(*),out)           - temporary array for material properties
c       khygro      (int,sc,out)             - hygro loading key
c       stress      (dp,ar(nTens),out)       - Cauchy stress at t+Dt
c       consJac     (dp,ar(nTens,nTens),out) - consistent Jacobian
c       strain      (dp,ar(nTens),out)       - total strain at t+Dt
c       siglv       (dp,ar(4),out)           - hygro load vector
c       siglvi      (dp,ar(4),out)           - initial stress load vector
c       eEnerPt     (dp,ar(7),out)           - energy at this point
c       nlkword     (int*8,ar(2,i),out)      - nonlinear material index
c             nlkword(1,1)                   - virtual starting position of plasticity
c             nlkword(2,1)                   - plasticity option
c             nlkword(1,2)                   - virtual starting position of the data
c             nlkword(2,2)                   - creep option
c
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
#if defined(MOIST)
#include "cont_info.inc"
      double precision eh
#endif
c
      integer          matRecIndex(*), matMWIndex(*), matWorkIndex(*),
     &                 nTens, elemId, kDomIntPt, kSectPt, kMatPoint,
     &                 matId, IntId, keycut, ireal
      integer*8        nlkword(W_nMatCol,W_nMatRow), 
     &                 keyplL, keycrL, keycpL, keyupL
      double precision matRec(*), matMW(*), matWork(*), prop(13),
     &                 strainInc(*),stress(*), consJac(nTens,*),
     &                 strain(*), tMatProp(*), siglv(4),
     &                 eEnerPt(7), siglvi(*), sdv(5),
     &                 coords(3), basis(3,3)
      double precision mref, moist, dmoist
  

      integer          khygro
      logical          iniFlag

      double precision ZERO
      parameter       (ZERO=0.0d0)


#if defined(DEBUG)
      call SubBeg('MaterialDriver')
#endif
      
      dmoist = ZERO
      
      keyplL = nlkword(1,1)
      keycrL = nlkword(1,2)
      keycpL = nlkword(1,3)
      keyupL = nlkword(1,4)
      
      sdv(1) = A_StraMax
      sdv(2) = A_MaxEqStra
      sdv(3) = A_CrMax

c --- get material index
      eEnerPt(1) = ZERO
      eEnerPt(4) = ZERO
      eEnerPt(5) = ZERO

c --- plasticity debug
c --- hygro strain calculations
#if defined(MOIST)
      call CalHygStrain(elemId, matId, 1, nTens,
     &                  moist, mref, matMW(1), matMWIndex(1), khygro)
      
      if (K_nMast.gt.0 .and. abs(eh).gt.TINY) then
          matMW(matMWIndex(W_iStrainHyg)) = 
     &                     matMW(matMWIndex(W_iStrainHyg)) + eh
          khygro = 1
      endif

#else      
      khygro = 0
      call v_zero(matMW(matMWIndex(W_iStrainHyg)),nTens)
#endif 

c --- solution variables update
      call SoluRetrieve (IntId, nTens, khygro,
     &                   nlkword(1,1), strainInc(1), stress(1),
     &                   matMW (1), matMWIndex(1))
c --- small strain linear elasticity
      call ElasticDriver (nTens, elemId, matId, moist,
     &                    matMW(1), matMWIndex(1),
     &                    consJac(1,1), prop(1), tMatProp(1),
     &                    khygro, siglv(1), kMatPoint)

c --- get initial stress
      if (iniFlag) then
          call StresIni (elemId, matId, nTens,
     &                   kDomIntPt, kSectPt, kMatPoint,
     &                   coords(1), basis(1,1), siglvi(1))
          call v1pv2_v1 (matMW(matMWIndex(W_iStress)), siglvi(1), nTens)
      endif

      if (keycrL.eq.0 .and. keyplL.eq.0 .and. keycpL.eq.0) then
          call ElasticUpdate (nTens,matMW(1),matMWIndex(1),
     &                        stress(1),eEnerPt(1))
      elseif (keycpL.eq.0 .and. keyupL.eq.0) then
c --- viscoplasticity/creep material model
          call PlCrDriver (keycut, moist,
     &                     nTens, elemId, matId, prop(1),
     &                     nlkword(1,1), matMW(1), matMWIndex(1),
     &                     matWork(1), matWorkIndex(1),
     &                     stress(1), consJac(1,1),
     &                     eEnerPt(1), eEnerPt(4), eEnerPt(5), sdv(1))
      elseif (keyupL.ne.0) then
          call uMatDriver (keycut, kMatPoint,
     &                     moist, dmoist,
     &                     nTens, elemId, matId,
     &                     nlkword(1,1), matMW(1), matMWIndex(1),
     &                     stress(1), consJac(1,1),
     &                     eEnerPt(1), eEnerPt(4), sdv(1))
      else
          call PlChabDriver (keycut, kMatPoint,
     &                       moist, dmoist,
     &                       nTens, elemId, matId, prop(1),
     &                       nlkword(1,1), matMW(1), matMWIndex(1),
     &                       matWork(1), matWorkIndex(1),
     &                       stress(1), consJac(1,1),
     &                       eEnerPt(1), eEnerPt(4), sdv(1))
      endif

      if (A_StraMax.lt.sdv(1) .or.
     &    A_MaxEqStra.lt.sdv(2) .or.
     &    A_CrMax.lt.sdv(3)) then
c
          call lockSMP (W_LockMaxStra)
          if (A_StraMax.lt.sdv(1)) A_StraMax = sdv(1)
          if (A_MaxEqStra.lt.sdv(2)) then
              A_MaxEqStra = sdv(2)
c --- candidate for deletion              
              if (A_MaxEqStra.gt.0.2d0) A_DelRealBeam = ireal
          endif
          if (A_CrMax.lt.sdv(3)) A_CrMax = sdv(3)  
          call unlockSMP (W_LockMaxStra)
c          if (A_StraMax.gt.A_StraLimit .and. A_IterNum.gt.2
c     &        .and. A_BisecNum.lt.2) keycut = 1
c          if (A_CrMax .gt. 0.2d0) keycut = 1
      endif

c --- update elastic strain energy
      matMW(matMWIndex(W_iMisc)) = eEnerPt(1)

c --- compute total strain

      call CalTotStrain (nTens, nlkword(1,1), khygro,
     &                   matMW(1), matMWIndex(1),
     &                   matRec(1), matRecIndex(1),
     &                   strain(1))

c --- update the material point variables to t+Dt
c --- store initial stress if necessarly
      if (iniFlag .and. A_SubStepNum.eq.1) then
          call PutMatRec (W_iStrainIni,kMatPoint,nTens,siglvi(1))
      endif


#if defined(DEBUG)
      call SubEnd('MaterialDriver')
#endif

      return
      end
