      subroutine BeamWorkIndx (lenmx,numElemNodes, off_secDataI,
     &                         off_RstNodeList,off_nCrnr,
     &                         off_work,off_twork,off_nResults,
     &                         off_matStrain,off_matStress,
     &                         off_secDataR,off_secResults,
     &                         off_matRecord,off_pResults,
     &                         off_postdata,off_ple,
     &                         off_peeq,off_pwrk,
     &                         off_thStrSecPt,off_tSecPt,
     &                         off_thStrSNDP,off_thStrSecNd,
     &                         off_cre,off_creq,off_crwk,
     &                         lwork, lenSecDataR,
     &                         nElemIntWork, nElemDbpWork,
     &                         nElemMisWork)
c --- primary function: calculate the length of element working memory for beam element
c
c  input arguments:
c     lenmx        (int,sc,in)       - length of materail data in SVR file
c     numElemNodes (int,sc,in)       - number of element nodes, including internal
c
c  output arguments:
c     off_xxxxxxxxxx (ar,sc,out)     - index for different working arrays
c     lwork          (int,sc,out)    - length of working array for stress output
c     lenSecDataR    (int,sc,out)    - length of working array for section data
c     nElemIntWork   (int,sc,out)    - number of integers in element working array
c     nElemDbpWork   (int,sc,out)    - number of double precision data in element working array
c     nElemMisWork   (int,sc,out)    - number of miscellenous data in element working array
c
#include "impcom.inc"
#include "constants.inc"
#include "sec_info.inc"
c
      integer   lenmx, lwork, numElemNodes,
     &          lenSecDataR, nElemIntWork, nElemDbpWork,
     &          nElemMisWork, len_postdata

      integer
     &          off_secDataI, off_RstNodeList, off_nCrnr,
     &          off_work, off_twork, off_nResults,
     &          off_matStrain, off_matStress, off_secDataR,
     &          off_secResults, off_matRecord, off_pResults,
     &          off_postdata, off_ple, off_peeq, off_pwrk,
     &          off_thStrSecPt, off_tSecPt, off_thStrSNDP,
     &          off_thStrSecNd, off_cre, off_creq, off_crwk


#if defined(DEBUG)
      call SubBeg('BeamWorkIndx')
#endif

      lwork            = lenmx + 1024
      len_postdata     = lenmx
      off_secDataI     = 1
      off_RstNodeList  = off_secDataI   + B_IntDatLen + 1000
      off_nCrnr        = off_RstNodeList+ B_MaxSecNds
      nElemIntWork     = off_nCrnr      + B_MaxSecNds
      off_work         = 1
      lenSecDataR      = B_MaxSecNds*3  + B_MaxSecPts/4
      off_twork        = off_work       + lwork
      off_nResults     = off_twork      + numElemNodes*lenSecDataR
      off_matStrain    = off_nResults   + 6*B_MaxSecNds*2
      off_matStress    = off_matStrain  + 3*B_MaxSecPts*2
      off_secDataR     = off_matStress  + 3*B_MaxSecPts*2
      off_secResults   = off_secDataR   + lenSecDataR*8
      off_matRecord    = off_secResults + B_NumResPts*B_MaxSecPts*8
      off_pResults     = off_matRecord  + lenmx
      off_postdata     = off_pResults   + 9*B_MaxSecNds*3
      off_ple          = off_postdata   + len_postdata
      off_peeq         = off_ple        + 3*B_MaxSecPts*2
      off_pwrk         = off_peeq       + 3*B_MaxSecPts*2
      off_thStrSecPt   = off_pwrk       + 3*B_MaxSecPts*2
      off_tSecPt       = off_thStrSecPt + B_MaxSecPts*2
      off_thStrSNDP    = off_tSecPt     + B_MaxSecPts*2
      off_thStrSecNd   = off_thStrSNDP  + B_MaxSecNds*2
      off_cre          = off_thStrSecNd + B_MaxSecNds*3
      off_creq         = off_cre        + 3*B_MaxSecPts*2
      off_crwk         = off_creq       + 3*B_MaxSecPts*2
      nElemDbpWork     = off_crwk       + 3*B_MaxSecPts*2
      nElemMisWork     = 75 + 9*B_MaxSecNds
      nElemIntWork     = nElemIntWork + 1
      nElemDbpWork     = nElemDbpWork + 1


#if defined(DEBUG)
      call SubEnd('BeamWorkIndx')
#endif

      return
      end
