      subroutine DynRotUpd (dun,vel,acel)
c --- primary function: calculate velocity as for now
c
c  input arguments:
c    dun      (dp,ar(3),in)    - time step rotational increment DOF vector
c
c  output arguments:
c    vel      (dp,ar(3),inout) - 1st time derivative of ut (velocity)
c    acel     (dp,ar(3),inout) - 2nd time derivative of ut (acceleration)
c
#include "impcom.inc"
#include "solu_info.inc"
#include "dyn_info.inc"
c
      double precision dun(3),vel(3),acel(3),vnew,alpha,con
      equivalence (alpha,I_Consts(2))

      if (A_Trans.eq.0) then
c --- a simple finite difference for structural dofs
         vel(1) = dun(1)/A_TimeInc
         vel(2) = dun(2)/A_TimeInc
         vel(3) = dun(3)/A_TimeInc
      else
         if (alpha.le.1.0d-6) then
c --- 2nd order newmark
c --- explicit (vel and acel are approximations)
             vnew = dun(1)/A_TimeInc
             con = 0.5d0/(A_TimeInc + A_TimeIno)
             acel(1) = con*(vnew - vel(1))
             acel(2) = con*(vnew - vel(2))
             acel(3) = con*(vnew - vel(3))
             vel(1) = vnew
             vel(2) = vnew
             vel(3) = vnew
         endif
      endif

      return
      end
