      function FormOutput ()
c --- primary function: result file writing leg
c
c  output arguments:
c     FormOutput  (int,sc,out) - statement number return
c
#include "impcom.inc"
#include "constants.inc"
#include "RST_file.inc"
#include "dof_info.inc"
#include "solu_info.inc"
#include "in_mem.inc"
c
      integer   SaveLoadSet, GetErr
      logical   IfEqvDP

      integer*8   lochedL, locrfL
      integer   FormOutput, i, nrf


#if defined(DEBUG)
      call SubBeg('FormOutput')
#endif

      FormOutput = 0

c --- set up the output masking

      M_PtrElemOutMask = PTRFTNNULL

c --- create the output controls
      call OutMaskControl (IfEqvDP(A_TimeCur,A_TimeEnd),R_OutMask)

c --- write out initial result file data
      call HedToRST
      i = SaveLoadSet (R_OutMask,A_TimeCur,A_SubStepNum,A_LoadStepNum,
     &                 A_CumIterNum,M_Disp(1),M_ForcTot(1),
     &                 locrfL,lochedL)

c --- check for error
      if (i.lt.0) then
          FormOutput = 999
          goto 999
      endif

      call OutIni

      call MakeElemOutput (D_NumEl,M_ElemOrd(1),
     &                     M_SoluVects(1),M_ForcTot(1))

      call MemFree (M_PtrElemOutMask)
      
      call SaveReact (R_OutMask,M_ForcTot(1),locrfL,nrf)

      call Buf_Flush (R_BlockNum)

  999 continue

#if defined(DEBUG)
      call SubEnd('FormOutput')
#endif

      return
      end

