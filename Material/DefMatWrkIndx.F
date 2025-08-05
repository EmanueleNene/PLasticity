      subroutine DefMatWrkIndx (nlkword, nTens, firstIndex, matWrkIndex)
c --- primary function: define material record index
c
c  input arguments:
c      nlkword    (intL,ar(2,*),in) - nonlinear  option  array
c      nTens      (int,sc,in)       - degree of freedom associated with element type
c
c  output arguments:
c      firstIndex   (int,sc,out)    - pointer to material work array
c      matWrkIndex  (int,ar(*),out) - material work array indices
c
#include "impcom.inc"
#include "constants.inc"
c
      integer    matWrkIndex(W_LenMatWrkIndx)
      integer*8  nlkword(W_nMatCol,W_nMatRow)
      integer    firstIndex, nTens, GetNumChabProp

      integer    nTensSq
      integer*8  keyplL, keycrL, keycpL, keyupL
      integer    matPtEqSize, matPtMiscSize, MiscSize
      parameter (matPtEqSize=2, matPtMiscSize=2)


      nTensSq = nTens * nTens
c --- nonlinear option index
      keyplL = nlkword(1,1)
      keycrL = nlkword(1,2)
      keycpL = nlkword(1,3)
      keyupL = nlkword(1,4)
     
      firstIndex = 1
      matWrkIndex(W_iStrain) = firstIndex
      matWrkIndex(W_iStress) = firstIndex + nTens
      matWrkIndex(W_iStrainInc) = firstIndex + nTens*2
      matWrkIndex(W_iStrainHyg) = firstIndex + nTens*3
      matWrkIndex(W_iMatMatrix) = firstIndex + nTens*4
      firstIndex = matWrkIndex(W_iMatMatrix) + nTensSq
c --- plastic strain
      if (keyplL.ne.0) then
          matWrkIndex(W_iStrainPl) = firstIndex
          matWrkIndex(W_iEqvPl) = firstIndex + nTens
          firstIndex = matWrkIndex(W_iEqvPl) + matPtEqSize
      endif
c --- creep strain
      if (keycrL.ne.0) then
          matWrkIndex(W_iStrainCr) = firstIndex
          matWrkIndex(W_iEqvCr) = firstIndex + nTens
          firstIndex = matWrkIndex(W_iEqvCr) + matPtEqSize
      endif

c --- miscellaneous record
      MiscSize = matPtMiscSize
      matWrkIndex(W_iMisc) = firstIndex
      firstIndex = firstIndex + MiscSize
c --- cyclic plasticity
      if (keyupL.ne.0) then
          matWrkIndex (W_iStateVar) = firstIndex
          firstIndex = firstIndex + max(9,nlkword(2,4))
      elseif (keycpL.ne.0) then
          matWrkIndex (W_iStateVar) = firstIndex
          firstIndex = firstIndex + GetNumChabProp(keycpL)*nTens
      endif


      
      return
      end
