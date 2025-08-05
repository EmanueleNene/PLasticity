      subroutine SoluRetrieve (iIntPt, nTens, khygro,
     &                         nlkword, strainInc, stress,
     &                         matWrk, matWrkIndex)
c --- primary function: retrieve the solution variables, such as stress, strain
c                       and others from previous sub-step
c
c  input arguments:
c      nTens        (int,sc,in)       - number of tensor components
c      keycrL       (int,sc,in)       - 1 creep, 0 otherwise
c      strainInc    (dp,ar(nTens),in) - incremental strain
c      matWrkIndex  (int,ar(*),in)    - material work area index
c      khygro       (int,sc,in)       - hygro flag
c      iniflag      (int,sc,in)       - initial stress flag
c      matWrkIndex  (int,ar(*),in)    - material MW Index
c
c  output arguments:
c      matWrk     (dp,ar(*),inout)    - appropriate strains in work areas
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
c
      double precision ZERO
      parameter       (ZERO = 0.d0)

      integer          nTens, khygro, matWrkIndex(*)
      integer*8        nlkword(W_nMatCol,W_nMatRow)
      double precision matWrk(*),
     &                 stress(nTens), strainInc(nTens)

      integer          iIntPt

      integer          leq
      integer*8        keyplL, keycrL
      double precision work(6)
      integer          matPtEqSize, matPtMiscSize
      parameter       (matPtEqSize=2, matPtMiscSize=2)


      keyplL = nlkword(1,1)
      keycrL = nlkword(1,2)

      call RetMatRec (W_iStress, iIntPt, nTens, stress(1) )
      call RetMatRec (W_iStrain, iIntPt,
     &                nTens, matWrk(matWrkIndex(W_iStrain)))

c --- strain increment
      call v_move (strainInc(1),
     &             matWrk(matWrkIndex(W_iStrainInc)), nTens)

      if (keyplL.gt.0) then
          call RetMatRec (W_iStrainPl, iIntPt,
     &                    nTens, matWrk(matWrkIndex(W_iStrainPl)) )
          call v1mv2_v1 (matWrk(matWrkIndex(W_iStrain)),
     &                   matWrk(matWrkIndex(W_iStrainPl)), nTens)
          leq = matPtEqSize
          call RetMatRec (W_iEqvPl, iIntPt,
     &                    leq, matWrk(matWrkIndex(W_iEqvPl)))
      endif

      if (keycrL.gt.0) then
          call RetMatRec (W_iStrainCr, iIntPt,
     &                    nTens, matWrk(matWrkIndex(W_iStrainCr)))
          call v1mv2_v1 (matWrk(matWrkIndex(W_iStrain)),
     &                   matWrk(matWrkIndex(W_iStrainCr)), nTens)
          leq = matPtEqSize
          call RetMatRec (W_iEqvCr, iIntPt,
     &                    leq, matWrk(matWrkIndex(W_iEqvCr)))
      endif

      leq = matPtMiscSize 

      call RetMatRec (W_iMisc, iIntPt,
     &                leq, matWrk(matWrkIndex(W_iMisc)) )

      if (khygro.gt.0) then
          call RetMatRec (W_iStrainHyg, iIntPt,nTens, work(1))
          call v1mv2_v1 (matWrk(matWrkIndex(W_iStrain)),
     &                   work(1), nTens)

          work(1) = matWrk(matWrkIndex(W_iStrainHyg)) - work(1)

          matWrk(matWrkIndex(W_iStrainHyg)) = work(1)

          if (A_KeyFirstIter.ne.1) then
              call v1mv2_v1 (matWrk(matWrkIndex(W_iStrainInc)),
     &                       work(1), nTens)
          endif
      endif

      call v1pv2_v1 (matWrk(matWrkIndex(W_iStrain)),
     &               matWrk(matWrkIndex(W_iStrainInc)), nTens)

      return
      end