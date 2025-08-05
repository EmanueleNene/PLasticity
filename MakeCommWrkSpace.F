      subroutine MakeCommWrkSpace (elemID, elChar, nUsrDof,
     &                             nElemIntWork, nElemDbpWork)
c --- primary function: calculate the length of element working memory space
c
c  input arguments:
c     elemID     (int,sc,in)          - element number
c     elChar     (int,ar(IELCSZ),in)  - element characteristic
c     nUsrDof    (int,sc,in)          - matrix and lv size
c
c  output arguments:
c     nElemIntWork (int,sc,out)       - number of integers in the work memory
c     nElemDbpWork (int,sc,out)       - number of double precision variables in the work memory
c
#include "impcom.inc"
#include "constants.inc"
c

      integer          elChar(IELCSZ), elemID, nUsrDof, nElemDbpWork,
     &                 nElemIntWork

      integer          nElStbDbWork, elStbIdx(W_LenStabIndx)


#if defined(DEBUG)
      call SubBeg('MakeCommWrkSpace')
#endif

      nElemIntWork = 0
      nElemDbpWork = 0

c --- index for stabilization

      if (elChar(NMSTVR).gt.0) then
          call DefStbIndx (elemID, elChar(1), nUsrDof, elStbIdx(1),
     &                     nElStbDbWork)
          nElemDbpWork = nElemDbpWork + nElStbDbWork
      endif


#if defined(DEBUG)
      call SubEnd('MakeCommWrkSpace')
#endif

      return
      end
