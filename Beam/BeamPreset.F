      subroutine BeamPreset (NumEl)
c --- primary function: allocate space for contact
c
c  input arguments:
c     NumEl (int,sc,in)          - number of elements
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "beam_hyg.inc"
c
      integer NumEl
#if defined(MOIST)
      PTRFTN  MemAlloc
      integer elmget,etyget,elem,nbeams,nnod,
     &        elData(EL_DIM),nodes(W_NdMax),elChar(IELCSZ),
     &        ityp,itypold


#if defined(DEBUG)
      call SubBeg('BeamPreset')
#endif

      M_PtrBeams = PTRFTNNULL
      M_PtrHygStrains = PTRFTNNULL
      H_NumBeams = 0
      M_PtrBeams = MemAlloc(NumEl,W_LenInt,'M_PtrBeams      ')
      call vi_zero (M_Beams(1),NumEl)
      
      nbeams = 0
c --- count beam elements
      itypold = 0
      do elem = 1,NumEl
         nnod = elmget(elem,elData(1),nodes(1))
         if (nnod.le.0) cycle
         if (elData(EL_DEAD).eq.1) cycle
         ityp = elData(EL_TYPE)
         if (ityp .ne. itypold) then
             itypold = ityp
             ityp = etyget(ityp,elChar(1))
         endif
c --- only beams are processed
         if (elChar(JBEAM).gt.0) then
             nbeams = nbeams + 1
             M_Beams(elem) = nbeams
         endif
      enddo

      if (nbeams.gt.0) then
          M_PtrHygStrains = MemAlloc(H_Col*H_Row*nbeams,
     &                               W_LenDouble,'M_PtrHygStrains ')
          H_NumBeams = nbeams
      else
          call erhandler ('BeamPreset',5000,2,
     &                    'No beam elements found',0.0d0,' ')          
      endif

#if defined(DEBUG)
      call SubEnd('BeamPreset')
#endif
#endif defined(MOIST)
      return
      end
