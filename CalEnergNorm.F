      subroutine CalEnergNorm (Fapp,FnrIn,uinc,lsval,lspar)
c --- primary function: compute the energy norm, du*(F-Fnr)
c
c  input arguments:
c     Fapp     (dp,ar(k7npc),in)   - applied load
c     FnrIn    (dp,ar(k7npc),in)   - newton-raphson restoring force
c     uinc     (dp,ar(k7npc),in)   - scaled displacement increment
c     lspar    (dp,sc,in)          - line search parameter that uinc has already
c                                    been scaled
c
c  output arguments:
c     lsval   (int,sc,out)         - computed e-norm
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "dof_info.inc"
c
      integer  i,j,k
      double precision  Fapp(*),FnrIn(*),uinc(*),lsval,lspar,df
      logical uconstraint

#if defined(DEBUG)
      call SubBeg('CalEnergNorm')
#endif

c --- compute dot product
      lsval = 0.0d0
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
                    df = Fapp(k) - FnrIn(k)
                    lsval = lsval + uinc(k)*df
                endif
            endif
         enddo

      enddo

c --- modify the value by the line search parameter
c --- this obtains the energy in terms of the unscaled displacement increment
      if (lspar.ne.0.0d0) lsval = lsval/lspar


#if defined(DEBUG)
      call SubEnd('CalEnergNorm')
#endif

      return
      end
