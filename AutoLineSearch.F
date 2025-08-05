      subroutine AutoLineSearch
c --- primary function: determine if line search should be on or off
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
#include "conv_stat.inc"
#include "in_mem.inc"
c
      integer i, ieqitm1
      double precision wnorm

#if defined(DEBUG)
      call SubBeg('AutoLineSearch')
#endif

c --- initialization
      ieqitm1 = A_IterNum - 1
      if (ieqitm1.eq.0) goto 999
      if (ieqitm1.eq.1) then
          do i = 1,5
             C_SavedInc(i) = 0.0d0
          enddo
      endif

c --- once line search is on, we let it on for the rest of the substep
      if (A_KeyLineSrch.eq.1) goto 999

c --- more than 3 
      if (ieqitm1.ge.3 .or. A_BisecNum.gt.1) then
	    A_KeyLineSrch = 1
          goto 999
      endif

      call CalEnergNormLS (M_ForcTot(1),M_ForcNR(1),M_IncI(1),wnorm)

      if (ieqitm1.ge.1.and.ieqitm1.le.4) then
          C_SavedInc(ieqitm1) = wnorm
      elseif (ieqitm1.ge.5) then
          C_SavedInc(5) = wnorm
      endif

c --- update work increments
      if (ieqitm1.ge.5) then
          do i = 1,4
             C_SavedInc(i) = C_SavedInc(i+1)
          enddo
      endif

c --- for 3 iterations the C_SavedInc oscillates, we turn
c --- line search on when it is far from convergence
c --- for contact situations:
      if (ieqitm1.ge.3) then
          if (C_SavedInc(1) .gt. C_SavedInc(2)) then
              if (C_SavedInc(2).lt.C_SavedInc(3)) A_KeyLineSrch = 1
          else
              if (C_SavedInc(2).gt.C_SavedInc(3)) A_KeyLineSrch = 1
          endif
      elseif (ieqitm1.eq.2) then
          if (C_SavedInc(1) .lt. C_SavedInc(2)) A_KeyLineSrch = 1
      endif

c --- set A_ConvLineSrch logical parameter for the solver
 999  if (A_KeyLineSrch.eq.1) then
          A_ConvLineSrch = .false.
      else
          A_ConvLineSrch = .true.
      endif
      if (A_KeyImpoDisp.eq.0) A_ConvLineSrch = .true.

c --- it is time to switch to full linesearch
      if (ieqitm1.ge.4 .or. A_BisecNum.gt.0) then
          A_KeyAutoLineSrch = 2
          A_KeyLineSrch = 1
      endif

#if defined(DEBUG)
      call SubEnd('AutoLineSearch')
#endif

      return
      end
