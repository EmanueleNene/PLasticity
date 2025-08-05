      subroutine RetMatRec0 (matRecordIndex, kMatRecPt, nTens, vect)
c --- primary function: retrieve initial material record data
c
c  input arguments:
c     matRecordIndex(in,sc,in)        - material record index
c     kMatRecPt     (in,sc,in)        - integration point where data
c                                     - to be inquired
c     nTens         (in,sc,in)        - number of components to be inquired
c     numWordsPt    (in,ar(*),in)     - number of material records
c
c  output arguments:
c     vect          (dp,ar(*),ou)     - inquired data
c
#include "impcom.inc"
#include "constants.inc"
#include "elem_wrk.inc"
#include "in_mem.inc"
c
      integer          iProcSMP
      integer          matRecordIndex, kMatRecPt, nTens
      double precision vect(*)
      integer          numWordsPt, i, k1, kk, iptr, dptr


      iptr = 1 + E_LenMatRecIndx*iProcSMP()

      if (M_MatRecI(iptr+matRecordIndex).gt.0) then
          numWordsPt = M_MatRecI(iptr)
          dptr = E_PtrMatRecDp0+ E_LenMatRec*iProcSMP()
          dptr = dptr - 1
          kk = dptr + (kMatRecPt-1) * numWordsPt
          k1 = kk + M_MatRecI(iptr+matRecordIndex) - 1

      else
          call v_zero(vect(1),nTens)
          goto 999
      endif

      do i = 1, nTens
         vect(i) = M_MatRecDp(k1+i)
      enddo


  999 continue

      return
      end

