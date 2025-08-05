      subroutine StabDriver (elemId, elChar, numDomainIntPts,
     &                       numElemNodes, nDofNode,
     &                       Gauss2, Gauss3, xBeg,
     &                       Gaussw2, Gaussw3, hbasis,
     &                       sectSummary, pIzz, pIyy, pIyz, pArea,
     &                       nGenStr, offset, lump,
     &                       nUsrDof, work, elDisp, kelout, eStiff,
     &                       fInt, elEnergy, ElRec, numSaveVar,
     &                       secStiff, gradMat)
c --- primary function: form mass matrix
c
c  input arguments:
c     elemId      (int,sc,in)          - element number
c     elChar      (int,ar(*),in)       - element characteristics
c     numDomainIntPts   (int,sc,in)    - number of domain integration points
c     xBeg        (dp,ar(nr,5),in)     - nodal coords (orig) and rot angles
c     numElemNodes(int,sc,in)          - number of element nodes
c     nUsrDof     (int,sc,in)          - matrix and LV size
c     Gauss2      (dp,ar(2),in)        - location of gaussian integration points
c     Gauss3      (dp,ar(3),in)        - location of gaussian integration points
c     Gaussw2     (dp,ar(2),in)        - weight factors of at gaussian points
c     Gaussw3     (dp,ar(3),in)        - weight factors of at gaussian points
c     hbasis      (dp,ar(3,3,3),in)    - base vectors at gaussian points
c     sectSummary (int,ar(20,8),in)    - section data
c     pIzz        (int,sc,in)          - Izz location in sectSummary
c     pIyy        (int,sc,in)          - Iyy location in sectSummary
c     pIyz        (int,sc,in)          - Iyz location in sectSummary
c     pArea       (int,sc,in)          - Area location in sectSummary
c     offset      (dp,ar(3,3,3),in)    - off set
c     lump        (int,sc,in)          - key for lumped mass
c     nDofNode    (int,sc,in)          - number of DOFs per node
c     elDisp      (dp,ar(nr,5),in)     - element nodal solution values
c
c  output arguments:
c     eStiff      (dp,ar(nr,nr),inout) - K matrix (A_MatrReq(1))
c     fInt        (dp,ar(nr),out)      - NR restoring F vector (A_MatrReq(6))
c     elEnergy    (dp,ar(10),out)      - element energies
c 
c
#include "impcom.inc"
#include "constants.inc"
#include "elem_wrk.inc"
#include "in_mem.inc"
#include "solu_info.inc"
c
      integer  iProcSMP

      integer
     &     elemId,numDomainIntPts, numElemNodes, pIzz, pIyy, pIyz,
     &     pArea, nGenStr, lump, nDofNode, nUsrDof,
     &     elChar(*), kelout(10), numSaveVar

      double precision
     &     Gauss2(2), Gauss3(3), xBeg(6,4), Gaussw2(2),Gaussw3(3),
     &     hbasis(3,3,3), sectSummary(20,8), cDensity, density,
     &     offset(3,3,3), work(*), elDisp(nUsrDof,*),
     &     eStiff(*), fInt(*), elEnergy(10), ElRec(numSaveVar),
     &     secStiff(nGenStr,nGenStr), gradMat(nGenStr,nDofNode*3)

      integer          nElStbDbWork, elStbIdx(W_LenStabIndx),
     &                 thisProc, idxStbFor, idxStbVar,
     &                 idxRecStbVar, nStbVar

c
      pointer          (pdeMass, eMass)
      pointer          (pdStbDpVar, StbDpVar)
      
      double precision eMass(nUsrDof,*), StbDpVar(*), dperr(3)

      integer          n, i, j, istar, jstar

#if defined(DEBUG)
      call SubBeg('MassForStab')
#endif
      
      if ((A_StabFirstStep.gt.0  .and.
     &     A_StabFirstStep.lt.3) .and. A_KeyMatUpd.eq.0) goto 999

c --- initialization
      call DefStbIndx (elemId, elChar(1), nUsrDof, elStbIdx(1),
     &                  nElStbDbWork)
      if (nElStbDbWork.gt.E_NumCommDp) then
          dperr(1) = elemId
          dperr(2) = nElStbDbWork
          dperr(3) = E_NumCommI
          call erhandler('MassForStab', 5000, 3, 'Element %I needs
     &                    dp work space length = %I, for stabilization,
     &                    which is larger than allocated length = %I.',
     &                    dperr(1),' ')
         goto 999
      endif
      thisProc = iProcSMP()
      pdeMass = M_PtrElemWrkDp + (E_NumWrkDp+elStbIdx(W_iStabMat)
     &                         - 1 + E_NumTotDp*thisProc)*W_AddrStrideDp

      call v_zero (eMass(1,1), nUsrDof*nUsrDof)

      cDensity = 1.0d0
      density = 1.0d0

      call BeamMass (numDomainIntPts, numElemNodes, nDofNode,
     &               Gauss2(1), Gauss3(1), xBeg(1,1), Gaussw2(1),
     &               Gaussw3(1), hbasis(1,1,1), sectSummary(1,1),
     &               pIzz, pIyy, pIyz, pArea, nGenStr, cDensity,
     &               density, offset(1,1,1), lump, 1,
     &               nUsrDof, secStiff(1,1), gradMat(1,1),
     &               work(1), eMass(1,1))

c --- neglect any translation and rotation coupling
      
      if (lump.ne.1) then
          do n = 1,numElemNodes
             istar = (n-1)*nDofNode
             jstar = 0
             do i = 1,numElemNodes
                do j = 1,3
                   eMass(istar+1, jstar+j+3) = 0.0d0
                   eMass(istar+2, jstar+j+3) = 0.0d0
                   eMass(istar+3, jstar+j+3) = 0.0d0
                   eMass(istar+4, jstar+j) = 0.0d0
                   eMass(istar+5, jstar+j) = 0.0d0
                   eMass(istar+6, jstar+j) = 0.0d0
                enddo
                jstar = jstar + nDofNode
             enddo
          enddo
      endif

      pdStbDpVar = M_PtrElemWrkDp + (E_NumWrkDp
     &                            +  E_NumTotDp*thisProc)*W_AddrStrideDp

      idxStbFor = elStbIdx(W_iStabForc)
      idxStbVar = elStbIdx(W_iStabVar)
      nStbVar = elChar(NMSTVR)
      idxRecStbVar = numSaveVar-nStbVar+1 

      call v_move (ElRec(idxRecStbVar), StbDpVar(idxStbVar), nStbVar)

c --- form the stabilization matrix
      call StabForm (elemId,nUsrDof, kelout(1),elDisp(1,A_RowIncS),
     &               eMass(1,1), elEnergy(1),
     &               StbDpVar(idxStbVar), StbDpVar(idxStbVar+1),
     &               StbDpVar(idxStbVar+2), eStiff(1), fInt(1),
     &               StbDpVar(idxStbFor))

c --- save stabilization variables
      if (A_KeyMatUpd.eq.1) then
          call v_move (StbDpVar(idxStbVar),ElRec(idxRecStbVar),nStbVar)
      elseif (A_StabFirstStep.eq.3 .and. A_IterNum.eq.1) then
          ElRec(idxRecStbVar) = StbDpVar(idxStbVar)
      endif

999   continue

#if defined(DEBUG)
      call SubEnd('MassForStab')
#endif

      return
      end
