      function UpdDisp (AllowBisec, dincmx, dofimx, kcut)
c --- primary function: update the displacements
c
c  input arguments:
c     AllowBisec   (int,sc,in)    - key indicates bisection is possible
c
c  output arguments:
c     kcut         (int,sc,out)    - kcut for bisection
c     dincmx       (dp,sc,out)     - max displacement increment
c     dofimx       (int,sc,out)    - dof at max displacement increment
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
#include "in_mem.inc"
#include "conv_stat.inc"
c
      integer            GetErr
      integer            UpdDisp, AllowBisec, kcut, dofimx, nondefnorm
      double precision   dincmx


#if defined(DEBUG)
      call SubBeg('UpdDisp')
#endif

      UpdDisp = 0
      nondefnorm = 0
c --- update the displacements
      call NR_Upd (M_IncI(1),M_Disp(1),M_IncS(1),M_Velo(1),M_Acce(1),
     &             dofimx,A_MaxDisp,dincmx,A_MaxDofInc)

c --- do not check if near convergence failure
      call CheckConvDisp (M_IncI(1),M_Disp(1),C_Disp,nondefnorm)
      
c --- dirty fix to skip the contact removal in UpdContConn
      if (nondefnorm.ne.0) C_nDeBo = -1

      if (abs(A_MaxDisp).gt.A_DispLimit) then
c --- DOF limit exceeded
          call PutErr (1)
      endif

      if (GetErr().eq.1 .or. nondefnorm.ne.0) then
c --- element formulation problem: bisection or exit
          if (AllowBisec.eq.1 .and. A_KeyFirstPass.eq.0) then
c --- try bisection
              call IfBisect (4,kcut)
              call PutErr (0)
              UpdDisp = 150
          else
c --- get out of here
              call StopAll (2)
              UpdDisp = 999
          endif
      endif


#if defined(DEBUG)
      call SubEnd('UpdDisp')
#endif

      return
      end
