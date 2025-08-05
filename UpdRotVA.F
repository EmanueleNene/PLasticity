      subroutine UpdRotVA (ut,dun,vel,acc,BackList)
c --- primary function: update rotational velocity and acceleration
c
c  input arguments:
c    ut       (dp,sc,in)       - total dof vector
c    dun      (dp,sc,in)       - time step increment dof vector
c
c  output arguments:
c    vel      (dp,sc,inout)    - 1st time derivative of ut (velocity)
c    acc      (dp,sc,inout)    - 2nd time derivative of ut (acceleration)
c
#include "impcom.inc"
#include "dof_info.inc"
#include "dyn_info.inc"
c
      double precision ut(*), dun(*), vel(*), acc(*)

      integer          BackList(D_BacLen), k, i, NodeExt, RotDofKey, 
     &                 iStart, iEnd
      double precision phi(3), omega(3), domega(3), R(3,3)

#if defined(DEBUG)
      call SubBeg('UpdRotVA')
#endif

c --- initial settings
      RotDofKey = 4
      k = RotDofKey
      call RangeSMP (D_BacLen,iStart,iEnd)
c --- node loop
      do i = iStart,iEnd
         NodeExt = BackList(i)
c --- dof loop
         if (ut(k).ne.HUGE.and.ut(k+1).ne.HUGE.and.ut(k+2).ne.HUGE) then
             phi(1) = dun(k)
             phi(2) = dun(k+1)
             phi(3) = dun(k+2)
             call GetRotMat (phi(1), R(1,1))
             omega(1) = -I_Prm(33)*vel(k) - I_Prm(32)*acc(k)
             omega(2) = -I_Prm(33)*vel(k+1) - I_Prm(32)*acc(k+1)
             omega(3) = -I_Prm(33)*vel(k+2) - I_Prm(32)*acc(k+2)
             domega(1) = -I_Prm(20)*vel(k) - I_Prm(21)*acc(k)
             domega(2) = -I_Prm(20)*vel(k+1) - I_Prm(21)*acc(k+1)
             domega(3) = -I_Prm(20)*vel(k+2) - I_Prm(21)*acc(k+2)
             vel(k) = I_Prm(31)*dun(k) +
     &                R(1,1)*omega(1) + 
     &                R(1,2)*omega(2) + 
     &                R(1,3)*omega(3)
             vel(k+1) = I_Prm(31)*dun(k+1) + 
     &                R(2,1)*omega(1) + 
     &                R(2,2)*omega(2) + 
     &                R(2,3)*omega(3)
             vel(k+2) = I_Prm(31)*dun(k+2) + 
     &                R(3,1)*omega(1) + 
     &                R(3,2)*omega(2) + 
     &                R(3,3)*omega(3)
             acc(k) = I_Prm(19)*dun(k) + 
     &                R(1,1)*domega(1) + 
     &                R(1,2)*domega(2) + 
     &                R(1,3)*domega(3)
             acc(k+1) = I_Prm(19)*dun(k+1) + 
     &                R(2,1)*domega(1) + 
     &                R(2,2)*domega(2) + 
     &                R(2,3)*domega(3)
             acc(k+2) = I_Prm(19)*dun(k+2) + 
     &                R(3,1)*domega(1) + 
     &                R(3,2)*domega(2) + 
     &                R(3,3)*domega(3)
         endif
      enddo

#if defined(DEBUG)
      call SubEnd('UpdRotVA')
#endif

      return
      end
