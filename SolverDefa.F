      subroutine SolverDefa
c --- primary function: initialize solution variables
c
#include "impcom.inc"
#include "solu_info.inc"
#include "conv_stat.inc"
c
      integer iwork
#if defined(DEBUG)
      call SubBeg('SolverDefa')
#endif
c --- set the stress calculation key for analyses
      if (A_KeyAutoStep.eq.1) then
          A_KeyBisect = 1
c --- set default nsub specifications here
          if (A_IniNSub.eq.0) then
c --- use A_IniNSub to define substep
              call erhandler('SolverDefa',5000,2,
     &          'The initial load step was chosen automatically',
     &           0.0d0,' ')
              A_IniNSub = 20
              if (A_MinNSub.eq.0) A_MinNSub = 4
              if (A_MaxNSub.eq.0) A_MaxNSub = 200000
          endif
c --- assign to all
          C_NSubIni = A_IniNSub
          C_NSubMax = A_MaxNSub
          C_NSubMin = A_MinNSub

          if (C_NSubMin.gt.C_NSubMax) then
              call erhandler('SolverDefa',5000,2,
     &          'Maximum number of substeps is lower than minimum.
     &          Values were swapped automatically.',
     &          0.0d0,' ')
              iwork = C_NSubMin
              C_NSubMin = C_NSubMax
              C_NSubMax = iwork
          endif


          if (C_NSubIni.lt.C_NSubMin) then
              call erhandler('SolverDefa',5000,2,
     &          'Initial number of substeps is lower than minimum.
     &           It was reset to minimum.',
     &           0.0d0,' ')
              C_NSubIni = C_NSubMin
          endif

          if (C_NSubIni.gt.C_NSubMax) then
              call erhandler('SolverDefa',5000,2,
     &          'Initial number of substeps is greater than maximum.
     &           It was reset to maximum.',
     &           0.0d0,' ')
              C_NSubIni = C_NSubMax
          endif

      endif

c --- initialize the 1st pass information
      A_KeyOut = 0
      A_CumIterNum = 0
      A_KeyPlast = 1
      A_KeyCreep = 2
      C_nAugm = 0

c --- initialize time variables
      A_TimeCur = 0.0d0
      A_TimeInc = 0.0d0
      A_TimeIno = 0.0d0

c --- initialize convergence parameters
      call v_zero (A_ConvCurVal(1),3)
      call v_zero (A_ConvPrvVal(1),3)

c --- other defaults
      if (A_SpecVer.eq.1) then
          A_MassLump  = 1
      else
          A_MassLump  = 0
      endif
      A_MaxLineSrchIter = 2
      A_LineSrchTol  = 0.50d0
      A_MinLineSrchVal = 0.05d0
      A_DispLimit = 1.0d6
      A_StraLimit = 0.03d0
      A_StraTol = 0.10d0
      A_CrLimit = 0.10d0
      A_CrTol = 0.10d0

#if defined(DEBUG)
      call SubEnd('SolverDefa')
#endif

      return
      end
