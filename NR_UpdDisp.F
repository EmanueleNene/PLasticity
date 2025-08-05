      subroutine NR_UpdDisp (dui,ut,dun,vel,acel,
     &                       glocmax,gumax,gdumax,gdnmax)
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
c    glocmax  (int,sc,out)     - DOF corresponding to dumax
c    gumax    (dp,sc,out)      - maximum DOF value
c    gdumax   (dp,sc,out)      - maximum incremental DOF value
c    gdumax   (dp,sc,out)      - maximum step-wise DOF value
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
#include "dof_info.inc"
#include "in_mem.inc"
#include "dyn_info.inc"
c
      integer  locmax,glocmax,i,j,k,RotDofKey,
     &         iStart,iEnd
      double precision
     &         ut(*),dun(*),dui(*),vel(*),acel(*),umax,dumax,dnmax,
     &         vnew,alpha,gumax,gdumax,gdnmax
      logical  mktiny(3)
      double precision ZERO, HALF
      parameter (ZERO = 0.0d0, HALF = 0.5d0)
      equivalence (alpha,I_Consts(2))

#if defined(DEBUG)
      call SubBeg('NR_UpdDisp')
#endif

c --- initial settings
      RotDofKey = 4
      call RangeSMP(D_BacLen,iStart,iEnd)
      k = D_NumDof*(iStart-1)
      dnmax = ZERO
      umax = ZERO
      dumax = ZERO
      
      do i = iStart,iEnd
c --- DOF loop
         j = 0
 100     j = j + 1
         if (j.le.D_NumDof) then

c --- increment the global DOF counter
            k = k + 1

c --- check for rotational update
            if (j.eq.RotDofKey) then
c --- all DOFs must be active
               if (dui(k).ne.HUGE .and. dui(k+1).ne.HUGE .and.
     &             dui(k+2).ne.HUGE) then
c --- check for inactive DOF becoming active
                  if (ut(k  ).eq.HUGE) ut(k  ) = ZERO
                  if (ut(k+1).eq.HUGE) ut(k+1) = ZERO
                  if (ut(k+2).eq.HUGE) ut(k+2) = ZERO

c --- check for zero pivots on each rotation
                  mktiny(1) = .false.
                  if (dui(k  ).eq.TINY) then
                      mktiny(1) = .true.
                      dui(k  ) = ZERO
                  endif
                  mktiny(2) = .false.
                  if (dui(k+1).eq.TINY) then
                      mktiny(2) = .true.
                      dui(k+1) = ZERO
                  endif
                  mktiny(3) = .false.
                  if (dui(k+2).eq.TINY) then
                      mktiny(3) = .true.
                      dui(k+2) = ZERO
                  endif

c --- now update the rotations
c --- total (incremental sum form):
                  ut(k  ) = ut(k  ) + dui(k  )
                  ut(k+1) = ut(k+1) + dui(k+1)
                  ut(k+2) = ut(k+2) + dui(k+2)
c --- total (pseudovector form):
                  call UpdRot (M_ERot(1,i),dui(k))
c --- increment over substep:
                  call UpdRot (dun(k),dui(k))
c --- transient updates
                  call DynRotUpd (dun(k),vel(k),acel(k))

                  if (abs(dui(k  )).gt.abs(dumax)) then
                      dumax = dui(k  )
                      locmax = k
                  endif
                  if (abs(dui(k+1)).gt.abs(dumax)) then
                      dumax = dui(k+1)
                      locmax = k+1
                  endif
                  if (abs(dui(k+2)).gt.abs(dumax)) then
                      dumax = dui(k+2)
                      locmax = k+2
                  endif
c --- get time step increment DOF value
                  if (abs(dun(k  )).gt.abs(dnmax))
     &                dnmax = dun(k  )
                  if (abs(dun(k+1)).gt.abs(dnmax))
     &                dnmax = dun(k+1)
                  if (abs(dun(k+2)).gt.abs(dnmax))
     &                dnmax = dun(k+2)
                  
                  if (abs(ut(k  )).gt.abs(umax)) umax = ut(k  )
                  if (abs(ut(k+1)).gt.abs(umax)) umax = ut(k+1)
                  if (abs(ut(k+2)).gt.abs(umax)) umax = ut(k+2)

c --- reset zero pivots
                  if (mktiny(1)) dui(k) = TINY
                  if (mktiny(2)) dui(k+1) = TINY
                  if (mktiny(3)) dui(k+2) = TINY
c --- now we've updated all 3 rotations
                  j = j + 2
                  k = k + 2
                  goto 100
               endif
            endif

c --- check for inactive DOF - skip it
            if (dui(k).eq.HUGE) then
                dun(k) = ZERO
                ut(k) = ZERO
                vel(k) = ZERO 
                if (M_PtrAcce.ne.PTRFTNNULL .and. A_NumSolVects.eq.5) 
     &              acel(k) = ZERO
                goto 100
            endif

c --- check for inactive DOF becoming active
            if (ut(k).eq.HUGE) ut(k) = ZERO
c --- check for "zero pivot" DOF
            mktiny(1) = .false.
            if (dui(k).eq.TINY) then
                mktiny(1) = .true.
                dui(k) = ZERO
            endif

c --- update total solution vector
            ut(k) = ut(k) + dui(k)
c --- update substep incremental solution vector
            dun(k) = dun(k) + dui(k)

c --- transient updates
            if (A_Trans.gt.0) then
c --- 2nd order newmark - structural
               if (alpha.gt.1.0d-6) then
c --- implicit
                   if (A_DynReverse.eq.0) then
                       call DynUpdFwd (dui(k),vel(k),acel(k),
     &                                 A_DynFirstPass)
                   else
                       call DynUpdBack (dui(k),vel(k),acel(k))
                   endif
               else
c --- explicit (vel and acel are approximations)
                   vnew = dun(k)/A_TimeInc
                   acel(k) = HALF*
     &                      (vnew - vel(k))/(A_TimeInc + A_TimeIno)
                   vel(k) = vnew
               endif
            else
c --- use a simple finite difference for structural DOFs
               vel(k) = dun(k)/A_TimeInc
            endif
c --- get time step increment DOF value
            if (abs(dun(k)).gt.abs(dnmax)) dnmax = dun(k)
            if (abs(ut(k)).gt.abs(umax)) umax = ut(k)
c --- get maximum values
            if (abs(dui(k)).gt.abs(dumax)) then
                dumax = dui(k)
                locmax = k
            endif
c --- reset "zero pivot" DOF
            if (mktiny(1)) dui(k) = TINY
c --- next node DOF
            goto 100
         endif
      enddo

c --- to avoid to many locks, we doublecheck
      if (abs(dumax).gt.abs(gdumax) .or. 
     &    abs(dumax).gt.abs(gdumax) .or.
     &    abs(umax).gt.abs(gumax)) then
          call lockSMP (W_LockElemInq)
          if (abs(dumax).gt.abs(gdumax)) gdumax = dumax
          if (abs(dnmax).gt.abs(gdnmax)) gdnmax = dnmax
          if (abs(umax).gt.abs(gumax)) then
              gumax = umax
              glocmax = locmax
          endif
          call unlockSMP (W_LockElemInq)
      endif    

      
#if defined(DEBUG)
      call SubEnd('NR_UpdDisp')
#endif

      return
      end