      subroutine CalBeamWrkSpace (elemID,elemChar,elData,
     &                            nElemIntWork,nElemDbpWork)
c --- primary function: calculate the length of element working memory
c
c  input arguments:
c     elemID   (int,sc,in)          - element number
c     elemChar (int,ar(IELCSZ),in)  - array of element type characteristics
c     elData   (int,ar(10),in)      - element attributes
c
c  output arguments:
c     nElemIntWork   (int,sc,out)   - number of integers in element working
c     nElemDbpWork   (int,sc,out)   - number of double precision data in element working array
c
#include "impcom.inc"
#include "constants.inc"
c
      integer
     &     elemID,elemChar(*),elData(*),nElemIntWork,
     &     nElemDbpWork, mRecIntPidx(W_LenMatRecIndx)

      integer
     &     off_secDataI, off_RstNodeList, off_nCrnr,
     &     off_work, off_twork, off_nResults,
     &     off_matStrain, off_matStress, off_secDataR,
     &     off_secResults, off_matRecord,
     &     off_pResults, off_postdata,
     &     off_ple, off_peeq, off_pwrk,
     &     off_thStrSecPt, off_tSecPt, off_thStrSNDP,
     &     off_thStrSecNd, off_cre, off_creq, off_crwk

      integer
     &     nmatMisc,numElemNodes,matId,lenmx,numWordsPt,
     &     lwork,lenSecDataR,nElemMisWork

      integer*8  nlkword(W_nMatCol,W_nMatRow)

#if defined(DEBUG)
      call SubBeg('CalBeamWrkSpace')
#endif

      numElemNodes = elemChar(NMNDNE)
      matId = elData(EL_MAT)

c --- material nonlinear properties index
      call GetMatKeys (matId,nlkword(1,1))
      call vi_zero (mRecIntPidx,W_LenMatRecIndx)
      call CalMatRecLen (nlkword(1,1), elemChar(1), elData(1),
     &                   lenmx, numWordsPt, nmatMisc, mRecIntPidx)

      call BeamWorkIndx (lenmx,numElemNodes, off_secDataI,
     &                   off_RstNodeList,off_nCrnr,
     &                   off_work,off_twork,off_nResults,
     &                   off_matStrain,off_matStress,
     &                   off_secDataR,off_secResults,
     &                   off_matRecord,off_pResults,
     &                   off_postdata,off_ple,
     &                   off_peeq,off_pwrk,
     &                   off_thStrSecPt,off_tSecPt,
     &                   off_thStrSNDP,off_thStrSecNd,
     &                   off_cre,off_creq,off_crwk,
     &                   lwork, lenSecDataR,
     &                   nElemIntWork, nElemDbpWork,
     &                   nElemMisWork)
      nElemDbpWork = nElemDbpWork + nElemMisWork

#if defined(DEBUG)
      call SubEnd('CalBeamWrkSpace')
#endif

      return
      end
