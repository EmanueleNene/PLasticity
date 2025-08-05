      subroutine CalEnergNormLS (Fapp,FnrIn,uinc,wnorm)
c --- primary function: compute the energy norm, du^t*(F-Fnr)
c                       for auto determination of line search on/off
c
c  input arguments:
c     Fapp     (dp,ar(k7npc),in)   - applied load
c     FnrIn    (dp,ar(k7npc),in)   - newton-raphson restoring force
c     uinc     (dp,ar(k7npc),in)   - scaled displacement increment
c
c  output arguments:
c     wnorm   (int,sc,out)         - computed e-norm
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "dof_info.inc"
c
      integer  i,j,k
      double precision  Fapp(*),FnrIn(*),uinc(*),wnorm,c2
      logical uconstraint

#if defined(DEBUG)
      call SubBeg('CalEnergNormLS')
#endif

c --- compute dot product
      wnorm = 0.0d0
      k = 0

c --- node loop
      do i = 1,D_BacLen
c --- DOF loop
         do j = 1,D_NumDof
c --- increment the DOF counter
            k = k + 1
c --- skip constrained DOFs
            uconstraint = btest(M_DofBits(k),W_DispBit)
            if (Fapp(k).eq.-HUGE) cycle
            if (Fapp(k).eq.HUGE)  uconstraint = .true.
            if (.not. uconstraint) then
c --- update energy norm
                if (uinc(k).ne.HUGE) then
c --- just use residual norm for now
                    c2 = Fapp(k) - FnrIn(k)
                    wnorm = wnorm + c2*c2
                endif
            endif
         enddo
      enddo
      wnorm = sqrt (wnorm)

#if defined(DEBUG)
      call SubEnd('CalEnergNormLS')
#endif

      return
      end
