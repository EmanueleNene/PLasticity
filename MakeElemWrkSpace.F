      subroutine MakeElemWrkSpace
c --- primary function: to create working in-memory space
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "elem_wrk.inc"
#include "dof_info.inc"
#include "solu_info.inc"
c
      integer  InqSMP, elmget, etyget

      PTRFTN   MemAlloc

      integer  itypeold, nElemIntWork, nElemDbpWork, nproc,
     &         ielem, elemID, numNodes, itype, nUsrDof, k, jtyp,
     &         elData(EL_DIM), nodes(W_NdMax), elChar(IELCSZ), length,
     &         nElemIntWorkCmn,nElemDbpWorkCmn

#if defined(DEBUG)
      call SubBeg('MakeElemWrkSpace')
#endif

c --- initialization
      itypeold = 0
      E_NumWrkI = 0
      E_NumWrkDp = 0
      nElemIntWork = 0
      nElemDbpWork = 0
      E_NumCommI = 0
      E_NumCommDp = 0
      nElemIntWorkCmn = 0
      nElemDbpWorkCmn = 0

      nproc = InqSMP(0)

c --- find the maximum length

      do 100 ielem = 1,D_NumEl
         elemID = M_ElemOrd(ielem)
         numNodes = elmget(elemID, elData(1), nodes(1))
         if (numNodes .le. 0) goto 100
         itype = elData(EL_TYPE)
         if (itype .ne. itypeold) then
             k = etyget(itype,elChar(1))
             k = k
             itypeold = itype
         endif
         nUsrDof = elChar(NMDFPN)*elChar(NMNDST)

c --- work space for each individual element
         jtyp = elChar(JETYP)
         if (jtyp.eq.W_eBeam2nd .or. jtyp.eq.W_eBeam3nd) then
             call CalBeamWrkSpace (elemID, elChar(1), elData(1),
     &                             nElemIntWork,nElemDbpWork)
         else
             goto 100
         endif

         if (E_NumWrkI.lt.nElemIntWork)  E_NumWrkI = nElemIntWork
         if (E_NumWrkDp.lt.nElemDbpWork) E_NumWrkDp = nElemDbpWork

c --- work space for material features of the element

         if (A_KeyStab.gt.0 .and. elChar(NMSTVR).gt.0) then
             call MakeCommWrkSpace (elemID, elChar(1),
     &                              nUsrDof, nElemIntWorkCmn,
     &                              nElemDbpWorkCmn)
            if (E_NumCommI.lt.nElemIntWorkCmn)
     &          E_NumCommI  = nElemIntWorkCmn
            if (E_NumCommDp.lt.nElemDbpWorkCmn)
     &          E_NumCommDp = nElemDbpWorkCmn
         endif

 100  continue

c --- create the integer index array
      E_NumTotI = E_NumWrkI+E_NumCommI
      if (E_NumTotI.gt.0) then
          length = E_NumTotI*nproc
          M_PtrElemWrkI = MemAlloc (length, W_LenInt,
     &                              'M_PtrElemWrkI   ')
          call vi_zero (M_ElemWrkI(1), length)
      endif


c --- create the double precision array
      E_NumTotDp = E_NumWrkDp+E_NumCommDp
      if (E_NumTotDp.gt.0) then
          length = E_NumTotDp*nproc
          M_PtrElemWrkDp = MemAlloc (length, W_LenDouble,
     &                               'M_PtrElemWrkDp  ')
          call v_zero (M_ElemWrkDp(1), length)
      endif


#if defined(DEBUG)
      call SubEnd('MakeElemWrkSpace')
#endif

      return
      end
