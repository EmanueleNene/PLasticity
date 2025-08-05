      function LineSearchDrv (CountLS,AllowBisec,dincmx,sdncmx,kcut)
c --- primary function: perform a line search
c
c  input arguments:
c     CountLS       (int,sc,in)   - line search counter
c     AllowBisec    (int,sc,in)   - key indicates bisection is possible
c     dincmx        (dp,sc,in)    - maximum displacement increment
c
c  output arguments:
c     sdncmx        (dp,sc,out)   - line search parameter
c     kcut          (int,sc,out)  - bisection parameter
c     LineSearchDrv (int,sc,out)  - return statement number
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "dof_info.inc"
#include "conv_stat.inc"
#include "solu_info.inc"
c
      integer   GetErr

      integer   LineSearchDrv, CountLS, AllowBisec, 
     &          lsitr, lspch, i, kcut

      double precision  dincmx, sdncmx, lspar, lsparo, lsval(3),
     &                  lsparr(3), enrm0, curBracket, oldBracket

#if defined(DEBUG)
      call SubBeg('LineSearchDrv')
#endif

      LineSearchDrv = 0
      do i = 1,3
         lsparr(i) = 0.0d0
         lsval(i) = 0.0d0
      enddo
      curBracket = 1.0d0
      oldBracket = 1.0d0

      lsitr = 1
      CountLS = CountLS + 1
c --- compute the initial energy norm (using the old fnr)
      lsparr(1) = 0.0d0
      lspar = 0.d0

      call CalEnergNorm (M_ForcTot(1),M_ForcNR(1),M_IncI(1),
     &                   enrm0,lspar)
      lsval(1) = enrm0
c --- now compute for the full increment
      lsparo = lspar
      lspar = 1.0d0
      lsparr(2) = 1.0d0

c >>> START LOOP
c --- initialize for load vector formation
  320 call LineSearchIni
c --- full DOF increment was applied below - readjust it
      if (lsitr.gt.1) then
          call AddIncDisp (lspar,lsparo,sdncmx)
      else
          sdncmx = dincmx
      endif

c --- form element load vectors
c --- only form the load vector
      A_InLineSrch = 1
      call FormMatr (D_NumEl,M_ElemOrd(1),M_SoluVects(1),.false.)
c --- add nodal forces to rhs for direct assembly
      if (D_NumEl.gt.0) call SumForce (M_ForcTot(1))

      A_InLineSrch = 0
      if (GetErr().eq.1 .or. C_CutMat.gt.0 .or. C_CutEle.gt.0) then
c --- error in element formulation: bisection or exit
          if (AllowBisec.eq.1 .and. A_KeyFirstPass.eq.0) then
c --- try bisection
              call IfBisect (5,kcut)
              call PutErr (0)
              LineSearchDrv = 150
              goto 999
          else
c --- print out error and get out of here
              if (GetErr().eq.0 .and. C_CutMat.eq.0) then
	            call DistortMessage(3)
              endif
              call StopAll (0)
              LineSearchDrv = 999
              goto 999
          endif
      endif

c --- form global load vectors
      A_KeyReform = 0

      if (GetErr().eq.1) goto 999

c --- calculate residual
      call CalEnergNorm (M_ForcTot(1),M_ForcNR(1),M_IncI(1),
     &                   lsval(2),lspar)

      lsparo = lspar

c --- decide if we need to continue and update parameter
      call CalLineSearchParam (lsitr,lsparr(1),lsval(1),enrm0,
     &                         lspar,curBracket,oldBracket,lspch)

      lsitr = lsitr + 1
      CountLS = CountLS + 1
      if (lspch.eq.1) goto 320

c <<< end LOOP

c --- set parameter to one if already converged
      if (C_Forc.eq.1 .and. C_Disp.eq.1) lspar = 1.0d0

c --- if parameter is one, set flag (cannot converge a substep
c --- with imposed displacement unless one equilibrium iter has lspar=1)
      if (lspar.gt.0.99999d0) A_ConvLineSrch = .true.
c --- adjust dof increment
      call AddIncDisp (lspar,lsparo,sdncmx)
      if (GetErr().eq.1) then
c --- dof limit exceeded: A_KeyBisect or exit
          if (AllowBisec.eq.1 .and. A_KeyFirstPass.eq.0) then
c --- try bisection
              call IfBisect (4,kcut)
              call PutErr (0)
              LineSearchDrv = 150
              goto 999
          else
c --- get out of here
              call StopAll (2)
              LineSearchDrv = 999
              goto 999
          endif
      endif

  999 continue

#if defined(DEBUG)
      call SubEnd('LineSearchDrv')
#endif

      return
      end
