      subroutine Predictor (dofimx,dincmx)
c --- primary function: perform displacement prediction for a substep
c
#include "impcom.inc"
#include "in_mem.inc"
#include "solu_info.inc"
#include "dof_info.inc"
c
      double precision  dincmx
      integer dofimx,i,node,idf,rotloc
      double precision  predfct
      logical IfEqvDP

#if defined(DEBUG)
      call SubBeg('Predictor')
#endif

      if (A_SubStepNum.gt.1 .and. A_Predict.gt.0 
     &    .and. A_BisecNum.eq.0) then
c --- find the start of ROTX
          rotloc = 4
c          if (A_Trans .eq. 0) then
c --- static/steady state
             if (IfEqvDP(A_TimeInc,A_TimeIno)) then
c --- time increments the same
                predfct = 1.0d0
                call v_move (M_IncS(1),M_IncI(1),D_ULen)
             else
c --- factor the increment by the ratio of the time steps
                if (A_TimeIno.eq.0.0d0) then
                    predfct = 1.0d0
                else
                    predfct = A_TimeInc/A_TimeIno
                endif
                predfct = min(1.5d0,predfct)
                call cv1_v2 (M_IncS(1),M_IncI(1),D_ULen,predfct)
             endif 
c          else
c --- dynamic/transient : use newmark predictor
c             call TransPred (M_IncI(1),M_Velo(1),M_Acce(1))
c          endif
c --- initialize the disp increment over the step
          call v_zero (M_IncS(1),D_ULen)
c --- NEW: start with zero prediction for rotation increment
          do i = 1,D_ULen
             call InqDof (i,node,idf)
             if (idf.ge.rotloc) M_IncI(i) = 0.0d0
          enddo
c --- update the displacements for the predictor
          call NR_Upd (M_IncI(1),M_Disp(1),M_IncS(1),M_Velo(1),
     &                 M_Acce(1),dofimx,A_MaxDisp,dincmx,A_MaxDofInc) 
c
      elseif (A_KeyFirstPass.eq.0) then
          call v_zero (M_IncI(1),D_ULen)
          call v_zero (M_IncS(1),D_ULen)
      endif

#if defined(DEBUG)
      call SubEnd('Predictor')
#endif
      return
      end     



      
      
     
      
      