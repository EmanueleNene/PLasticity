      subroutine AddIncDisp (lspar,lsparo,dincmx)
c --- primary function: add the incremental disp
c
c  input arguments:
c     D_ULen     (int,sc,in)         - total number of dof in this analysis
c     lspar    (dp,sc,in)          - current line search parameter
c     lsparo   (dp,sc,in)          - last line search parameter
c
c  internal variables:
c     dincmx   (dp,sc)             - maximum displacement increment
c     dofimx   (int,sc)            - dof where the maximum increment occured
c
#include "impcom.inc"
#include "constants.inc"
#include "dof_info.inc"
#include "in_mem.inc"
#include "solu_info.inc"
c
      integer  dofimx,i
      double precision dincmx,lspar,lsparo,con

#if defined(DEBUG)
      call SubBeg('AddIncDisp')
#endif

      if (lspar.ne.lsparo) then

c --- negate the increment in displacement for removal
         do i = 1,D_ULen
            if (M_IncI(i).ne.HUGE .and. M_IncI(i).ne.TINY) then
                M_IncI(i) = -M_IncI(i)
            endif
         enddo

c --- set flag for reversal of Newmark accels and velocities
         A_DynReverse = 0
         if (A_Trans.gt.0 .and. A_IterNum.eq.1 
     &       .and. A_Predict.eq.0) then
             A_DynReverse = 1
         endif
         
c --- pass the reverse increment to the geometric update
         call NR_Upd (M_IncI(1),M_Disp(1),M_IncS(1),M_Velo(1),M_Acce(1),
     &                dofimx,A_MaxDisp,dincmx,A_MaxDofInc)

c --- displacement from the last iter is recovered
c     now process uinc to add new fraction of last solution

         con = lspar/lsparo

         do i = 1,D_ULen
            if (M_IncI(i).ne.HUGE .and. M_IncI(i).ne.TINY) then
                M_IncI(i) = -con*M_IncI(i)
            endif
         enddo

c --- pass the new increment to the geometric update
         call NR_Upd (M_IncI(1),M_Disp(1),M_IncS(1),M_Velo(1),M_Acce(1),
     &                dofimx,A_MaxDisp,dincmx,A_MaxDofInc)
      endif
      
      A_DynReverse = 0

#if defined(DEBUG)
      call SubEnd('AddIncDisp')
#endif

      return
      end
