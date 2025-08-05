      subroutine IncoreMatIndx (mRecIntPidx, nRecIntPidx, numWordsPt)
c --- primary function: update material record index to in-core memory
c
c  input arguments:
c     mRecIntPidx  (int, ar(*),in)   - material record index
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "elem_wrk.inc"
c
      integer  numWordsPt, nRecIntPidx, mRecIntPidx(nRecIntPidx)
      integer  iProcSMP, ptri, ptr, i, i1, i2

      ptri = E_PtrMatRecI  + E_LenMatRec*iProcSMP()
      i2   = ptri

      if(nint(M_MatRecDp(i2)) .le. 0.0d0) then
         M_MatRecDp(i2) = E_MatRecAdd
         i2  = i2 + 1
         M_MatRecDp(i2) = numWordsPt
         do i = 1, nRecIntPidx
            i2 = i2 + 1
            M_MatRecDp(i2) = mRecIntPidx(i)
         enddo
      endif
      ptr = 1 + E_LenMatRecIndx*iProcSMP()
      i1 = ptr - 1
      i2 = ptri
      do i = 1, nRecIntPidx + 1
         i1 = i1 + 1
         i2 = i2 + 1
         M_MatRecI(i1) = nint(M_MatRecDp(i2))
      enddo
         
      return
      end
