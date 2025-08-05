      subroutine DynUpdFwd (dui,vel,acel,DynFirstPass)
c --- primary function: newmark time integration updating
c
c  input arguments:
c    dui      (dp,sc,in)       - displacement increment
c
c  output arguments:
c    vel      (dp,sc,out)      - velocity
c    acel     (dp,sc,out)      - acceleration
c
#include "impcom.inc"
#include "dyn_info.inc"
c
      integer DynFirstPass
      double precision dui,vel,acel,vold,aold

c --- second order newmark (hyperbolic) update
      vold = vel
      aold = acel
      if (DynFirstPass .eq. 1) then
c --- first pass
          acel = I_Prm(19)*dui - I_Prm(20)*vold - I_Prm(21)*aold
          vel  = vold + I_Prm(7)*aold + I_Prm(8)*acel
      else
c --- subsequent passes
          acel = aold + I_Prm(19)*dui
          vel  = vold + I_Prm(8)*(acel-aold)
      endif

      return
      end
