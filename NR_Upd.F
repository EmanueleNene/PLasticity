      subroutine NR_Upd (dui,ut,dun,vel,acel,locmax,umax,dumax,dnmax)
c --- primary function: update total solution vectors for newton raphson
c
c  input arguments:
c    dui      (dp,sc,in)       - iteration increment DOF vector
c
c  output arguments:
c    ut       (dp,sc,inout)    - total primary DOF vector
c    dun      (dp,sc,inout)    - time step increment DOF vector
c    vel      (dp,sc,inout)    - 1st time derivative of ut (velocity)
c    acel     (dp,sc,inout)    - 2nd time derivative of ut (acceleration)
c    locmax   (int,sc,out)     - DOF corresponding to dumax
c    umax     (dp,sc,out)      - maximum DOF value
c    dumax    (dp,sc,out)      - maximum incremental DOF value
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
#include "in_mem.inc"
#include "cont_info.inc"
c
      external NR_UpdDisp
      integer  locmax, nprocs, InqSMP
      double precision
     &         ut(*),dun(*),dui(*),vel(*),acel(*),umax,dumax,dnmax
      double precision ZERO
      parameter (ZERO = 0.0d0)

#if defined(DEBUG)
      call SubBeg('NR_Upd')
#endif

c --- initial settings
      if (K_nMast.gt.0.and.A_IterNum.gt.0) then
          call UpdPiloNd (K_nPerFounI,K_nFoun,K_NdCoord(1),
     &                    M_BackList(1),K_PiloNd(1),
     &                    K_FounI(1),M_SoluVects(1))
      endif

      nprocs = InqSMP(0)
      locmax = 0
      umax = ZERO
      dumax = ZERO
      dnmax = ZERO
      if (nprocs.gt.1) then
          call InfoSMP (1,0)
          call ResumeSMP
          call ForkSMP9 (NR_UpdDisp,dui(1),ut(1),dun(1),vel(1),acel(1),
     &                   locmax,umax,dumax,dnmax)
          call SuspendSMP
      else
          call NR_UpdDisp (dui(1),ut(1),dun(1),vel(1),acel(1),
     &                     locmax,umax,dumax,dnmax)
      endif

c --- update initial transient update key
      if (A_DynReverse.eq.0) then
          A_DynFirstPass = 0
      else
          A_DynFirstPass = 1
      endif

#if defined(DEBUG)
      call SubEnd('NR_Upd')
#endif

      return
      end