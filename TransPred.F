      subroutine TransPred (dui,vel,acel)
c --- primary function: newmark predictor logic
c
c  input arguments:
c    vel      (dp,sc,in)     - 1st time derivative of ut (velocity)
c    acel     (dp,sc,in)     - 2nd time derivative of ut (acceleration)
c
c  output arguments:
c    dui      (dp,sc,inout)  - iteration increment dof vector
c
#include "impcom.inc"
#include "constants.inc"
#include "dof_info.inc"
#include "in_mem.inc"
#include "solu_info.inc"
#include "dyn_info.inc"
c
      integer conchk
      integer i,j,k,node,kerr
      double precision dui(*),vel(*),acel(*),alpha,work(4)

      equivalence (alpha,I_Consts(2))

c --- initial settings
      k = 0
c --- node loop
      do 200 i = 1,D_BacLen
c --- get expanded node number
         node = M_BackList(i)
         if (node.le.0) goto 200
c --- dof loop
         do 100 j = 1,D_NumDof
c --- increment the dof counter
            k = k + 1
c --- check for inactive dof
            if (dui(k).eq.HUGE) goto 100
c --- check for "zero pivot" dof
            if (dui(k).eq.TINY)  goto 100
c --- transient updates
            if (alpha.gt.1.0d-6) then
c --- implicit only
                dui(k) = A_TimeInc*( 
     &                   vel(k) + A_TimeInc*(0.5d0-alpha)*acel(k))
            endif

c --- is this degree of freedom a specified displacement?
            if (.not.(btest(M_DofBits(k),W_DispBit))) goto 100
            if (.not.(btest(M_DofBits(k),W_TranBit))) goto 100
c --- adjust for imposed displacements
            if (conchk(node,j).ne.2) then
                call GetD (node,j,work(1))
                if (work(1).ne.0.0d0) 
     &              dui(k) = work(1) - M_Disp(k)
            endif
  100    continue
  200 continue

      return
      end
