      subroutine WorkSpacesIni
c --- primary function: create element work space
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "elem_wrk.inc"
c
      integer          InqSMP
      PTRFTN           MemAlloc

      integer          nTotalSpace_HEAP,
     &                 nTotalSpace_HEAPI,
     &                 nWrkSpace_MatRec0
      integer          nproc, ptrMatRecBeg

      PTRFTN           lastloc

#if defined(DEBUG)
      call SubBeg('WorkSpacesIni')
#endif

      nproc = InqSMP(0)

      E_MatRecAdd = W_LenMatRecIndx + 3

c --- material record work space

      call MakeMatWrkSpace (nWrkSpace_MatRec0)

      E_LenMatRec = nWrkSpace_MatRec0 + E_MatRecAdd + 1

c --- total work space in HEAP
      nTotalSpace_HEAP = E_LenMatRec*3 + 1
c --- multiple processors
      nTotalSpace_HEAP  = nTotalSpace_HEAP*nproc

      M_PtrMatRecDp = MemAlloc(nTotalSpace_HEAP,W_LenDouble,
     &                        'M_PtrMatRecDp   ')

      call v_zero (M_MatRecDp(1),nTotalSpace_HEAP)

c --- work space pointer
      lastloc = 1

c --- material record pointer
      ptrMatRecBeg   = lastloc
c
      E_PtrMatRecDp  = ptrMatRecBeg + E_MatRecAdd
      E_PtrMatRecDp0 = E_PtrMatRecDp + E_LenMatRec*nproc

      E_PtrMatRecI   = ptrMatRecBeg
      E_PtrMatRecI0  = E_PtrMatRecI + E_LenMatRec*nproc
c --- FYI
c     lastloc        = E_PtrMatRecI + E_LenMatRec*nproc*2

c --- element index HEAP work space
      E_LenMatRecIndx = W_LenMatRecIndx + 1

c --- total index HEAP work space
      nTotalSpace_HEAPI = E_LenMatRecIndx
      nTotalSpace_HEAPI = nTotalSpace_HEAPI*nproc

c --- create index HEAP work space
      M_PtrMatRecI = MemAlloc(nTotalSpace_HEAPI,W_LenInt,
     &                        'M_PtrMatRecI    ')
      call vi_zero (M_MatRecI(1),nTotalSpace_HEAPI)

c --- material record index
      lastloc = 1 + E_LenMatRecIndx*nproc
      
      
#if defined(DEBUG)
      call SubEnd('WorkSpacesIni')
#endif

      return
      end
