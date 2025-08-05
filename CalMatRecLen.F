      subroutine CalMatRecLen (nlkword, elChar, elData,
     &                         lenmx, numwd,
     &                         nmatMisc, mRecIntPidx)
c --- primary function: determine the maximum length of material record required
c
c  input arguments:
c     elChar        (int,ar(*),in)       - array of element type characteristics
c
c  output arguments:
c      nlkword      (intL,ar(2,*),inout)  - material flags
c      lenmx        (int,sc,out)          - numwd * number of integration points
c      numwd        (int,sc,out)          - total double precision words of
c                                           material record at integration point
c      nmatMisc     (int,sc,out)          - size of miscellaneous material record
c      mRecIntPidx  (int,ar(*),out)       - material record pointer at integration points
c
#include "impcom.inc"
#include "constants.inc"
c
      integer          elChar(*), mRecIntPidx(*),
     &                 lenmx, numwd, nmatMisc, elData(EL_DIM)
      integer*8        nlkword(W_nMatCol,W_nMatRow)

      integer          GetSecChars, numIntg

#if defined(DEBUG)
      call SubBeg('CalMatRecLen')
#endif

c --- initialize
      numwd = 0
      nmatMisc = 0
      lenmx = 0
            
      call DefMatRecIndx (elChar(1), nlkword, mRecIntPidx,
     &                    nmatMisc, numwd)

      numIntg = GetSecChars(elData(EL_REAL),1)
c --- was checked with IsQuadBeam before
      if (elChar(JETYP).eq.W_eBeam3nd) numIntg = 2*numIntg

      lenmx = numwd * numIntg

  999 continue
#if defined(DEBUG)
      call SubEnd('CalMatRecLen')
#endif

      return
      end
