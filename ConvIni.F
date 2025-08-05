      subroutine ConvIni
c --- primary function: convergence params initializations
c
#include "impcom.inc"
#include "solu_info.inc"
#include "constants.inc"
#include "conv_stat.inc"

      integer cnvfld(4), k, j
      integer kloc(4)
      data    kloc / 1, 3, 2, 4/

#if defined(DEBUG)
      call SubBeg('ConvIni')
#endif

      C_KeyMono = 0

c --- determine global options from the currently active DOF options      
      cnvfld (1) = 1         ! F
      cnvfld (3) = 1         ! M
      if (A_SpecVer.gt.0) then
          cnvfld (2) = 2     ! U
          if (A_SpecVer.eq.3 .and. C_nPivots.eq.0) cnvfld(2) = 0     ! U
          if (A_SpecVer.eq.4) cnvfld(2) = 0                          ! U
          cnvfld (4) = 0     ! R
      else
          cnvfld (2) = 0     ! U
          cnvfld (4) = 0     ! R
      endif
      
c --- l2-norm out-of-balance convergence checking
      do k = 1,4
         j = kloc(k) 
         A_ConvTol(j) = 0.005d0*dble(cnvfld(j))
         A_ConvNorm(j) = 2
         A_ConvRef(j) = 0.0d0
         A_ConvMinVal(j) = 0.0d0
      enddo
      
      if (C_nAugm.gt.0.and.C_IterCtrlG.gt.8.and.C_IterCtrlG.lt.12)
     &    C_IterCtrlG = 12

c --- for the first 2 substeps of a load step, we allow contact to iterate
c     until iter == A_IterMax to establish contact

      A_KeyLineSrch = 0
      if (A_KeyAutoLineSrch.eq.2) A_KeyLineSrch = 1
      A_InLineSrch = 0

c --- initialize convergence keys
      A_KeyConv = 0
      C_Forc = 0
      C_Disp = 0

c --- initialize contact status keys
      C_nPen = 0
      C_nPenIni = 0
      C_ConnUpd = 0
            
      C_ContChang = C_nConF + C_nDeBo + C_nConA

      C_nConF = 0            
      C_nDeBo = 0
      C_nConA = 0
      C_nFraBo = 0
      C_nPivots = 0
      
#if defined(DEBUG)
      call SubEnd('ConvIni')
#endif

      return
      end


