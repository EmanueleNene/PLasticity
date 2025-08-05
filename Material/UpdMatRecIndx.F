      subroutine UpdMatRecIndx (kMatPoint, nTens, khygro, nlkword,
     &                          strain, stress, matWrkIndex, matWrk)
c --- primary function: update the material record
c
c  input arguments:
c      nTens    (int,sc,in)          - number of tensor components
c      khygro   (int,sc,in)          - hygro flag
c      strain   (dp,ar(nTens),in)    - total strain
c      stress   (dp,ar(nTens),in)    - stress
c      matWrkIndex(int,ar(*),in)      - array indexing into matWrk
c      matWrk    (dp,ar(*),in)        - converged material variables
c      keypl    (int,sc,loc)         - >0 plasticity, 0 otherwise
c      keycr    (int,sc,loc)         - >0 creep     , 0 otherwise
c
#include "impcom.inc"
#include "constants.inc"
c
      integer          kMatPoint,
     &                 nTens, matWrkIndex(*),
     &                 khygro
      double precision
     &                 strain(nTens), stress(nTens),
     &                 matWrk(*)
c
      integer          lenwk
c
      double precision tmpwrk(6)
      integer          matPtEqSize, matPtMiscSize
      parameter       (matPtEqSize=2, matPtMiscSize=2)
      integer*8        nlkword(W_nMatCol,W_nMatRow)
      
      
      lenwk = matPtMiscSize
      call PutMatRec (W_iMisc, kMatPoint,
     &                lenwk, matWrk(matWrkIndex(W_iMisc)))
c
      if (khygro.gt.0) then
          lenwk = nTens
          call RetMatRec (W_iStrainHyg,kMatPoint,lenwk, tmpwrk(1))
          call v1pv2_v1 (tmpwrk(1),
     &                   matWrk (matWrkIndex (W_iStrainHyg)), lenwk)
          call PutMatRec (W_iStrainHyg, kMatPoint,lenwk, tmpwrk(1))
      endif

      lenwk  = nTens
      call PutMatRec (W_iStrain, kMatPoint, lenwk, strain(1))
      call PutMatRec (W_iStress, kMatPoint, lenwk, stress(1))

      if (matWrkIndex(W_iStrainPl).gt.0) then
          lenwk  = nTens
          call PutMatRec (W_iStrainPl, kMatPoint,
     &                    lenwk, matWrk(matWrkIndex(W_iStrainPl)))
c
          lenwk  = matPtEqSize
          call PutMatRec (W_iEqvPl, kMatPoint,
     &                    lenwk, matWrk(matWrkIndex(W_iEqvPl)))
      endif
c
      if (matWrkIndex(W_iStrainCr).gt.0) then
          lenwk  = nTens
          call PutMatRec (W_iStrainCr, kMatPoint,
     &                    lenwk, matWrk(matWrkIndex(W_iStrainCr)))
          lenwk  = matPtEqSize
          call PutMatRec (W_iEqvCr  , kMatPoint,
     &                    lenwk, matWrk(matWrkIndex(W_iEqvCr  )))
      endif

      return
      end
