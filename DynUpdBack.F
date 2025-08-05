      subroutine DynUpdBack (dui,vel,acel)
c -- primary function: reverse newmark update
c
c  input arguments:
c    dui      (dp,sc,in)       - displacement
c
c  output arguments:
c    vel      (dp,sc,inout)    - velocity
c    acel     (dp,sc,inout)    - acceleration
c
#include "impcom.inc"
#include "dyn_info.inc"
c
      double precision dui,vel,acel,acl1bk

c --- second order newmark (hyperbolic) update
      acl1bk = acel
      acel = (-I_Prm(19)*dui - I_Prm(20)*vel + 
     &       (I_Prm(20)*I_Prm(8)-1.0d0)*acel)/(I_Prm(21) - 
     &        I_Prm(20)*I_Prm(7))
      vel = vel - I_Prm(8)*acl1bk - I_Prm(7)*acel

      return
      end
