      subroutine SetDispPointers
c --- primary function: set up displacement pointers
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "dof_info.inc"
#include "solu_info.inc"
c
#if defined(DEBUG)
      call SubBeg('SetDispPointers')
#endif

c --- define the element displacement pointers from the global pointers
      A_RowDisp = 0
      A_RowIncI = 0
      A_RowIncS = 0
      A_RowVelo = 0
      A_RowAcce = 0
      A_NumSolVects = 0
      if (M_PtrDisp .ne. PTRFTNNULL) then
          A_RowDisp = (M_PtrDisp-M_PtrSoluVects)/(D_ULen*W_AddrStrideDp)
     &              + 1
          A_NumSolVects = A_NumSolVects + 1
      endif
      if (M_PtrIncI .ne. PTRFTNNULL) then
          A_RowIncI = (M_PtrIncI-M_PtrSoluVects)/(D_ULen*W_AddrStrideDp)
     &              + 1
          A_NumSolVects = A_NumSolVects + 1
      endif
      if (M_PtrIncS .ne. PTRFTNNULL) then
          A_RowIncS = (M_PtrIncS-M_PtrSoluVects)/(D_ULen*W_AddrStrideDp)
     &              + 1
          A_NumSolVects = A_NumSolVects + 1
      endif
      if (M_PtrVelo .ne. PTRFTNNULL) then
          A_RowVelo = (M_PtrVelo-M_PtrSoluVects)/(D_ULen*W_AddrStrideDp)
     &              + 1
          A_NumSolVects = A_NumSolVects + 1
      endif
      if (M_PtrAcce .ne. PTRFTNNULL) then
          A_RowAcce = (M_PtrAcce-M_PtrSoluVects)/(D_ULen*W_AddrStrideDp)
     &              + 1
          A_NumSolVects = A_NumSolVects + 1
      endif


#if defined(DEBUG)
      call SubEnd('SetDispPointers')
#endif

      return
      end
