      subroutine DefStbIndx (elId, elChar, nUsrDof, elStbIdx,
     &                        nElStbDbWork)
c --- primary function: build the kinematic work space index for stabilization
c
c  input arguments:
c     elId         (int,sc,in)           - element number
c     elChar       (int,ar(IELCSZ),in)   - element characteristic
c     nUsrDof      (int,sc,in)           - matrix and lv size
c
c  output arguments:
c     elStbIdx     (int,ar(*),out)       - index for stabilization
c     nElStbDbWork (int,sc,out)          - number of dp variables in the work memory for stabilization
c
#include "impcom.inc"
#include "constants.inc"
c
      integer          elId, elChar(IELCSZ), nUsrDof,
     &                 elStbIdx(W_LenStabIndx), nElStbDbWork

      integer          firstWord

#if defined(DEBUG)
      call SubBeg('DefStbIndx')
#endif

c --- initialization
      call vi_zero (elStbIdx(1), W_LenStabIndx)

      firstWord = 1
      elId = elId

      elStbIdx(W_iStabMat) = firstWord
      firstWord = firstWord + nUsrDof*nUsrDof
      elStbIdx(W_iStabForc) = firstWord
      firstWord = firstWord + nUsrDof
      elStbIdx(W_iStabVar) = firstWord
      firstWord = firstWord + elChar(NMSTVR)

      nElStbDbWork = firstWord - 1


#if defined(DEBUG)
      call SubEnd('DefStbIndx')
#endif

      return
      end


