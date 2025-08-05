      subroutine OutMaskControl (lTime,R_OutMask)
c --- primary function: store the stress output controls
c
#include "impcom.inc"
#include "constants.inc"
#include "dof_info.inc"
#include "in_mem.inc"
c
      PTRFTN    MemAlloc
      logical  lTime       ! if true, set up controls
c                          ! if false, set array to 0
      integer  R_OutMask   ! global output mask
      integer  leng

#if defined(DEBUG)
      call SubBeg('OutMaskControl')
#endif

c --- allocate space for stress output controls
      leng = max(D_NumEl,D_NumNd)
      M_PtrElemOutMask = MemAlloc(leng,W_LenInt,'M_PtrElemOutMask')

c --- initialize the output mask to the default values
      if (lTime) then
          call SetOutDefault (D_NumEl,D_NumNd,M_ElemOutMask(1),leng)
      else
          call vi_zero (M_ElemOutMask(1),leng)
          R_OutMask = 0
      endif

c --- set the output mask accordingly
      call FormOutMask (M_ElemOutMask(1),leng,D_NumEl,D_NumNd)

#if defined(DEBUG)
      call SubEnd('OutMaskControl')
#endif

      return
      end
