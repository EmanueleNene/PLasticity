      subroutine BeamElemDriver (elemID,elemChar,elData,eomask,eConn,
     &                           locsvrL,A_MatrReq,nUsrDof,xBeg,eleDisp,
     &                           kelout,eStiff,eMass,gStiff,fExt,fInt,
     &                           elemVolume,elemMass,elemCG,eEnergy,
     &                           edindxL,lcerstL)
c --- primary function: prepare input for beam driver
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "elem_wrk.inc"
c
      integer*8 locsvrL, edindxL(25), lcerstL
      integer   elemID,elemChar(*),elData(*),eomask,eConn(*),
     &          A_MatrReq(*),nUsrDof,kelout(*),
     &          numElemNodes,numSectPts, mRecIntPidx(W_LenMatRecIndx)
      integer*8 nlkword(W_nMatCol,W_nMatRow)

      double precision
     &     xBeg(*),eleDisp(*),elemCG(*),eEnergy(*),
     &     eStiff(*),eMass(*),gStiff(*),fExt(*),fInt(*),
     &     elemVolume,elemMass
      double precision
     &     eU(21),incDisp(21),secStiff(49),tGenStress(14),
     &     tGenTh(7),gradMat(150),dperr(3)
c
      integer matId
      integer
     &     off_secDataI, off_RstNodeList, off_nCrnr,
     &     off_work, off_twork, off_nResults,
     &     off_matStrain, off_matStress, off_secDataR,
     &     off_secResults, off_matRecord, off_pResults, off_postdata,
     &     off_ple, off_peeq, off_pwrk, off_thStrSecPt, off_tSecPt,
     &     off_thStrSNDP, off_thStrSecNd, off_cre, off_creq,
     &     off_crwk
      integer num_int_words, num_dp_words, num_ms_words,
     &        nDofNode,nGenStr, iProcSMP,
     &        lenmx, numWordsPt, lwork,ptr,lenSecDataR,
     &        nmatMisc, thisProc, GetSecChars

      pointer (pMI,MI)
      pointer (pMD,MD)
      pointer (pmisc,misc)
      integer MI(*)
      double precision MD(*), misc(*)

#if defined(DEBUG)
      call SubBeg('BeamElemDriver')
#endif

      nDofNode = 6
      nGenStr = 6

c --- account for internal nodes
      numElemNodes = elemChar(NMNDNE)

c --- material nonlinear properties index
      matId = elData(EL_MAT)
      call GetMatKeys (matId,nlkword(1,1))

      call vi_zero (mRecIntPidx(1),W_LenMatRecIndx)
      call CalMatRecLen (nlkword(1,1), elemChar(1), elData(1),
     &                   lenmx, numWordsPt,
     &                   nmatMisc, mRecIntPidx(1))
      
      call BeamWorkIndx (lenmx,numElemNodes, off_secDataI,
     &                   off_RstNodeList, off_nCrnr,
     &                   off_work, off_twork, off_nResults,
     &                   off_matStrain, off_matStress,
     &                   off_secDataR, off_secResults,
     &                   off_matRecord, off_pResults,
     &                   off_postdata, off_ple,
     &                   off_peeq, off_pwrk,
     &                   off_thStrSecPt, off_tSecPt,
     &                   off_thStrSNDP, off_thStrSecNd,
     &                   off_cre, off_creq, off_crwk,
     &                   lwork, lenSecDataR,
     &                   num_int_words, num_dp_words,
     &                   num_ms_words)

      if (num_dp_words+num_ms_words.gt.E_NumWrkDp) then
          dperr(1) = elemID
          dperr(2) = num_dp_words+num_ms_words
          dperr(3) = E_NumWrkDp
          call erhandler ('BeamElemDriver', 5000, 3, 'Element %I needs
     &                     double precision working space length = %I,
     &                     which is larger than created length = %I.',
     &                     dperr(1),' ')
          goto 999
      endif
      if (num_int_words.gt.E_NumWrkI) then
          dperr(1) = elemID
          dperr(2) = num_int_words
          dperr(3) = E_NumWrkI
          call erhandler ('BeamElemDriver', 5000, 3, 'Element %I needs
     &                     integer working space length = %I,
     &                     which is larger than created length = %I.',
     &                     dperr(1),' ')
          goto 999
      endif
      thisProc = iProcSMP()
      ptr = E_PtrMatRecDp + E_LenMatRec*thisProc
      pMD = M_PtrElemWrkDp + E_NumTotDp*thisProc*W_AddrStrideDp
      pMI = M_PtrElemWrkI + E_NumTotI*thisProc*W_AddrStrideI
      pmisc = pMD + num_dp_words*W_AddrStrideDp

      numSectPts = GetSecChars(elData(EL_REAL),1)

c --- initialization is no longer needed
      call v_zero(MD(1),num_dp_words)
      call vi_zero(MI(1),num_int_words)

      call BeamElement (elemID,elemChar,elData,eomask,eConn,locsvrL,
     &                  A_MatrReq,nUsrDof,xBeg,eleDisp,kelout,eStiff,
     &                  eMass,gStiff,fExt,fInt,elemVolume,elemMass,
     &                  elemCG,eEnergy,edindxL,lcerstL,numElemNodes,
     &                  lenmx, numWordsPt, nmatMisc,MI(off_secDataI),
     &                  MI(off_RstNodeList),MI(off_nCrnr),lwork,
     &                  MD(off_work),MD(off_twork),MD(off_nResults),
     &                  MD(off_matStrain),MD(off_matStress),
     &                  MD(off_secDataR),lenSecDataR,numSectPts,
     &                  MD(off_secResults),M_MatRecDp(ptr),
     &                  MD(off_pResults),MD(off_postdata),
     &                  MD(off_ple),MD(off_peeq),
     &                  MD(off_pwrk),MD(off_thStrSecPt),
     &                  MD(off_tSecPt),MD(off_thStrSNDP),
     &                  MD(off_thStrSecNd),nDofNode,nGenStr,
     &                  eU(1),incDisp(1),secStiff(1),
     &                  tGenTh(1),gradMat(1),tGenStress(1),
     &                  MD(off_cre),MD(off_creq),MD(off_crwk),
     &                  mRecIntPidx(1),nlkword(1,1),misc(1))

  999 continue

#if defined(DEBUG)
      call SubEnd('BeamElemDriver')
#endif

      return
      end
