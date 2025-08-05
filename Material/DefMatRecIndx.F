      subroutine DefMatRecIndx (elChar, nlkword, mRecIntPidx,
     &                          nmatMisc, numwd)
c --- primary function: define material record index at integration point
c
c  input arguments:
c     nlkword       (int*8,ar(2,*),in) - material flags
c
c  output arguments:
c     mRecIntPidx   (int,ar(*),out)    - material record pointer at
c                                        integration points
c     nmatMisc      (int,sc,out)       - size of miscellaneous material record
c
#include "impcom.inc"
#include "constants.inc"
c
      integer            numwd, nmatMisc, nTens, elChar(*),
     &                   mRecIntPidx (W_LenMatRecIndx)
      integer*8          nlkword(W_nMatCol,W_nMatRow)

      integer            firstIndex, GetNumChabProp
      integer*8          keyplL, keycrL, keycpL, keyupL

      integer    matPtEqSize, matPtMiscSize
      parameter (matPtEqSize=2,matPtMiscSize=2)

      nTens = elChar(NCOMPN)
      keyplL = nlkword(1,1)
      keycrL = nlkword(1,2)
      keycpL = nlkword(1,3)
      keyupL = nlkword(1,4)
c --- this will not work. we need this space to form the output
c      if ((keyplL+keycrL+keycpL).eq.0 .and. elChar(KYOP1).ne.1) then
c          nmatMisc = 0
c          numwd = 0
c          goto 999
c      endif
      firstIndex = 1
      mRecIntPidx (W_iStrain) = firstIndex
      mRecIntPidx (W_iStress) = firstIndex + nTens
      mRecIntPidx (W_iStrainHyg) = firstIndex + nTens*2
      firstIndex = firstIndex + nTens*3
      if (keyplL.ne.0) then
          mRecIntPidx (W_iStrainPl) = firstIndex
          mRecIntPidx (W_iEqvPl) = firstIndex + nTens
          firstIndex = firstIndex + nTens + matPtEqSize
      endif
      if (keycrL.ne.0) then
          mRecIntPidx (W_iStrainCr) = firstIndex
          mRecIntPidx (W_iEqvCr) = firstIndex + nTens
          firstIndex = firstIndex + nTens + matPtEqSize
      endif

c --- initial strain
      mRecIntPidx (W_iStrainIni) = firstIndex
      firstIndex = firstIndex + nTens

      nmatMisc = matPtMiscSize
c --- miscellaneous record
      mRecIntPidx (W_iMisc) = firstIndex
      firstIndex = firstIndex + nmatMisc
      
      if (keyupL.ne.0) then
          mRecIntPidx (W_iStateVar) = firstIndex
          firstIndex = firstIndex + max(9,nlkword(2,4))
      elseif (keycpL.ne.0) then
          mRecIntPidx (W_iStateVar) = firstIndex
          firstIndex = firstIndex + GetNumChabProp(keycpL)*nTens
      endif
      
      numwd = firstIndex - 1
999   continue

      return
      end