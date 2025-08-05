      subroutine UpdSvrRec (key, elemId, edindx, lenmx)
c --- primary function: update SVR records
c
c  input arguments:
c     key         (int,sc,in)         0 = get record
c                                     1 = conditionally put record
c                                    -1 = always put record
c     elemId      (int,sc,in)         - element number
c     lenmx       (int,sc,in)         - length of material record
c
c  output arguments:
c     edindx      (int,ar(20),out)    - element result data file indexes
c
#include "impcom.inc"
#include "in_mem.inc"
#include "solu_info.inc"
#include "constants.inc"
#include "elem_wrk.inc"
c
      integer          iProcSMP

      integer          elemId, key, lenmx, edindx(20)
      integer          matRecSet, nsvr, nsvr0
      parameter       (matRecSet = 2)

      integer          ptr0, ptri, i
      double precision dperr(1)

#if defined(DEBUG)
      call SubBeg('UpdSvrRec')
#endif

      ptr0 = E_PtrMatRecI0 + E_LenMatRec*iProcSMP()
      ptri = E_PtrMatRecI  + E_LenMatRec*iProcSMP()

      nsvr = lenmx + E_MatRecAdd

      if (key .eq. 0) then
          nsvr0 = nsvr
          call GetSVR (edindx(1), matRecSet, nsvr0,
     &                 M_MatRecDp(ptri))
         
          if (nsvr0 .lt. nsvr) then
              nsvr = nsvr0
              M_MatRecDp(ptri) = 0.0d0
          endif

          M_MatRecDp(ptri+E_MatRecAdd-1) = dble(nsvr)
         
          do i = 0,nsvr-1
             M_MatRecDp(ptr0+i) = M_MatRecDp(ptri+i)
          enddo

      elseif (key .eq. 1) then
          nsvr = nint(M_MatRecDp(ptri+E_MatRecAdd-1))
          if (nsvr .gt. lenmx) then
              ptri = E_PtrMatRecI  + E_LenMatRec*iProcSMP()
              ptr0 = E_PtrMatRecI0 + E_LenMatRec*iProcSMP()
          else
              ptri = E_PtrMatRecDp  + E_LenMatRec*iProcSMP()
              ptr0 = E_PtrMatRecDp0 + E_LenMatRec*iProcSMP()
          endif
          if (A_KeyFirstPass.eq.1) then
              call StoreSVR (edindx(1), matRecSet, nsvr ,
     &                       M_MatRecDp(ptr0))
          else
              call StoreSVR (edindx(1), matRecSet, nsvr ,
     &                       M_MatRecDp(ptri))
          endif
      elseif (key .eq. -1) then
          call StoreSVR (edindx(1), matRecSet, nsvr ,
     &                   M_MatRecDp(ptri))
      else
c --- should not be here, error
          dperr(1) = dble(elemId)
          call erhandler ('UpdSvrRec',5000,3,
     &        'Element %I is trying to access unknown SVR record',
     &         dperr(1),' ')
      endif

#if defined(DEBUG)
      call SubEnd('UpdSvrRec')
#endif

      return
      end
