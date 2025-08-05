      subroutine GetStabData  (elId, elChar, nUsrDof, key, nData,
     &                            sData)
c --- primary function: obtain sData related for stabilization
c
c  input arguments:
c     elId     (int,sc,in)           - element number
c     elChar   (int,ar(IELCSZ),in)   - element characteristic
c     nUsrDof  (int,sc,in)           - matrix and lv size
c     key      (int,sc,in)           - tell what sData is asked
c                                      = -1, all
c                                      =  1, stabilization matrix only
c                                      =  2, stabilization force only
c                                      =  3, energies only
c     nData    (int,sc,in)           - number of sData requested
c
c
c  output arguments:
c     sData     (dp,ar(*),out)       - returned sData
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "elem_wrk.inc"
c
      integer          iProcSMP

      integer          elChar(IELCSZ), elId, nUsrDof, key, nData

      double precision sData(*)

      integer          nElStbDbWork, elStbIdx(W_LenStabIndx),
     &                 idxStbMat, idxStbFor, idxStbVar, thisProc,
     &                 nUsrDof2, iData

      pointer         (pdStbDpVar, StbDpVar)

      double precision StbDpVar(*), dperr(3)

c --- initialization

      call v_zero (sData(1), nData)
      call DefStbIndx (elId, elChar(1), nUsrDof, elStbIdx(1),
     &                      nElStbDbWork)
      if (nElStbDbWork.gt.E_NumCommDp) then
          dperr(1) = elId
          dperr(2) = nElStbDbWork
          dperr(3) = E_NumCommI
          call erhandler ('GetStabData', 5000, 3, 'Element %I needs
     &                     double precision working space length = %I,
     &                     for stabilization, which is larger than
     &                     allocated length = %I.',dperr(1),' ')
          goto 999
      endif
      thisProc = iProcSMP()
      pdStbDpVar = M_PtrElemWrkDp + (E_NumWrkDp
     &           + E_NumTotDp*thisProc)*W_AddrStrideDp
      idxStbMat = elStbIdx(W_iStabMat)
      idxStbFor = elStbIdx(W_iStabForc)
      idxStbVar = elStbIdx(W_iStabVar)

      nUsrDof2 = nUsrDof*nUsrDof
      if (key.eq.1) then
          call v_move (StbDpVar(idxStbMat), sData(1), nUsrDof2)
          iData = nUsrDof2
      elseif (key.eq.2) then
          call v_move (StbDpVar(idxStbFor), sData(1), nUsrDof)
          iData = nUsrDof
      elseif (key.eq.3) then
          call v_move (StbDpVar(idxStbVar), sData(1), 3)
          iData = 3
      elseif (key.eq.-1) then
          call v_move (StbDpVar(idxStbMat), sData(idxStbMat), nUsrDof2)
          call v_move (StbDpVar(idxStbFor), sData(idxStbFor), nUsrDof)
          call v_move (StbDpVar(idxStbVar), sData(idxStbVar), 3)
          iData = nUsrDof2+nUsrDof+3
      endif

      if (iData.ne.nData) then
          dperr(1) = elId
          dperr(2) = nData
          dperr(3) = iData
          call erhandler ('GetStabData', 5000, 3, 'Element %I should
     &                     have %I variables for stabilization,
     &                     but only %I are available.',dperr(1),' ')
      endif

 999  continue

      return
      end


