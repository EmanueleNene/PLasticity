      subroutine PlChabDriver (keycut, kMatPoint, moist, dmoist,
     &                         nTens, elemId, matId, prop,
     &                         nlkword, matMW, matMWIndex,
     &                         matWork, matWorkIndex,
     &                         stress, dsdePl,
     &                         sedEl, sedPl, sdv)
c --- primary function: cyclic plasticity driver
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
c
      integer          nxiMax  
      parameter       (nxiMax = 5)

      integer          keycut, nTens, elemId, matId, GetNumChabProp,
     &                 matMWIndex(*), matWorkIndex(*), kMatPoint
     
      integer*8        nlkword(W_nMatCol,W_nMatRow), keycpL
      double precision matMW(*), matWork(*),
     &                 dsdePl(nTens,nTens), 
     &                 stress(nTens), sedEl, sedPl, sdv(*)
      double precision moist, dmoist, prop(13)
      double precision propcp(1+nxiMax*3), propcp0(1+nxiMax*3)

      integer          i1, i2, i3, i4, i5, i6, lenwk,
     &                 w1, w2, w3, w4, w5, w6, w7, w8, w9
      double precision time, dtime

      dtime = A_TimeInc
      time = A_TimeCur
      keycpL = nlkword(1,3)
      
     
c --- define indices
      i1 = matMWIndex(W_iStrain)
      i2 = matMWIndex(W_iStrainPl)
      i3 = matMWIndex(W_iEqvPl)
      i4 = matMWIndex(W_iStress)
      i5 = matMWIndex(W_iMisc)
      i6 = matMWIndex(W_iStateVar)
      
      w1 = matWorkIndex(W_WrkIndx1)
      w2 = matWorkIndex(W_WrkIndx2)
      w3 = matWorkIndex(W_WrkIndx3)
      w4 = matWorkIndex(W_WrkIndx4)
      w5 = matWorkIndex(W_WrkIndx5)
      w6 = matWorkIndex(W_WrkIndx6)
      w7 = matWorkIndex(W_WrkIndx7)
      w8 = matWorkIndex(W_WrkIndx8)
      w9 = matWorkIndex(W_WrkIndx9)
c --- store state variables (back stresses)
      if (keycpL.gt.0) then
          lenwk = GetNumChabProp(keycpL)*nTens
          call RetMatRec (W_iStateVar,kMatPoint,lenwk,matMW(i6))
      endif
      
      call ChabocheDriver (keycut,time, dtime,moist, dmoist,
     &                     nTens, elemId, matId, prop(1),
     &                     nxiMax, propcp0(1), propcp(1),
     &                     matMW(i1), matMW(i2),
     &                     matMW(i3), matMW(w8), matMW(i6),
     &                     matMW(i5),
     &                     matWork(w1), matWork(w2), matWork(w3),
     &                     matWork(w4), matWork(w5), matWork(w6),
     &                     matWork(w7), matWork(w9),
     &                     matMW(i4), stress(1), dsdePl(1,1) ,
     &                     nlkword(1,1), sedEl, sedPl,
     &                     sdv(1))

      if (keycpL.gt.0) then
          call PutMatRec (W_iStateVar,kMatPoint,lenwk,matMW(i6))
      endif

      
      return
      end
