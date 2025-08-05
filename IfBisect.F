      subroutine IfBisect (biloc,kcut)
c --- primary function: auto load bisection routine
c
c  input arguments:
c     biloc   (int,sc,in)      - key that indicates type of solver call
c
c  output arguments:
c     kcut    (int,sc,out)     - key that indicates that bisection will now be done
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "dof_info.inc"
#include "conv_stat.inc"
#include "solu_info.inc"
c
      logical  IfEqvDP, mintime

      integer  biloc, iott, kcut, biprnt, k, GetWriteUnit,
     &         newcutback, checkmono,  neqitr80percent,
     &         kStbCut
      double precision  con, dtnext, cut, stcrit, crcrit

#if defined(DEBUG)
      call SubBeg('IfBisect')
#endif

      mintime = IfEqvDP(A_TimeInc,A_MinTimeStep) .or. 
     &                 (A_TimeInc.lt.A_MinTimeStep)
      iott = GetWriteUnit()
      kcut = 0
      cut = 1.0d0
      newcutback = 0
      neqitr80percent = A_IterMax-2
      if (neqitr80percent.le.0) neqitr80percent = 1

c --- check different bisection cases
      if (biloc .eq. 1) then
c --- bisection: poor convergence prediction
c --- excluding the case that user does not want it,
c --- or the 1st 2 substep of contact establishment
c --- solu convergence pattern is not monotonic
         checkmono = 0
c         if (A_KeyAutoLineSrch.eq.0 .and. C_KeyMono.gt.1 
c     &       .and. A_SubStepNum.gt.1) then
c --- without linesearch - bisect upon first divergence
c             checkmono = 1
c       else
        if (C_KeyMono.ge.4 .and. A_SubStepNum.gt.1) then
             if (.not.mintime) then
                 checkmono = 1
             endif
         endif

         if (checkmono.eq.0 .and. C_nAugm.gt.0 
     &       .and. A_IterNum.lt.neqitr80percent) goto 788
     
         if (A_IterNum.ge.6) then
             k = A_IterMax
             k = k + 1
             if (C_nIterEsti.gt.k .and. C_nIterEsti1b.gt.k .and.
     &           C_nIterEsti2b.gt.k .and. A_IterNum.ge.C_IterCtrlR) then
                 biprnt = 1
                 kcut = 1
                 cut = 0.5d0
                 A_DynReverse = 0
             endif
         endif

         if (checkmono.gt.0.and.A_IterNum.ge.C_IterCtrl0) then
             biprnt = 1
             kcut = 1
             cut = 0.45d0
             newcutback = 1
             A_DynReverse = 0
         endif

  788    C_nIterEsti2b = C_nIterEsti1b
         C_nIterEsti1b = C_nIterEsti
      elseif (biloc .eq. 2) then
c --- bisection: negative pivot
         kcut = 1
         biprnt = 2
         if (A_Predict.gt.0) then
             cut = 1.0d0
             A_Predict = 0
         else
             cut = 0.5d0
         endif
         if (A_SubStepNum.eq.1) cut = 0.5d0
      elseif (biloc .eq. 3) then
c --- bisection: used all equilibrium iterations
         kcut = 1
         biprnt = 3
         cut = 0.5d0
      elseif (biloc .eq. 4) then
c --- bisection: dof limit exceeded
         newcutback = 1
c --- it is over predict displacement, turn off predictor try again
c --- no time cutback yet.
         kcut = 1
         biprnt = 4
         if (A_Predict.gt.0) then
             cut = 1.0d0
             A_Predict = 0
         else
             cut = 0.5d0
         endif
         if (A_SubStepNum.eq.1) cut = 0.5d0
      elseif (biloc .eq. 5) then
c --- bisection: error in element formulation
         kcut = 1
         biprnt = 5
         cut = 0.35d0
      elseif (biloc .eq. 7) then
c --- bisection: contact status change
         kcut = 1
         biprnt = 9
         cut = 0.25d0
      else
c --- bisection: plasticity or creep limits exceeded
         if (A_KeyPlast .gt. 0) then
             stcrit = A_StraLimit * (1.0d0 + A_StraTol)
             if (A_StraMax .gt. stcrit) then
                 kcut = 1
                 biprnt = 6
                 cut = A_StraLimit / A_StraMax
c --- if A_StraMax > A_StraLimit here I am going to use golden CUT=0.88d0
c     to address cases cut = 1.0d0 where no bisection would occur
                 if (cut .ge. 0.99d0) then
                     cut = 0.88d0
                 endif
             endif
         endif
c --- bisection:  creep ratio control
         if (A_KeyCreep.eq.1) then
c --- explicit creep
             crcrit = A_CrLimit * ( 1.d0 + A_CrTol )
             if (A_CrMax .gt. crcrit) then
                 con = A_CrLimit / A_CrMax
c --- same principle as plasticity
                 if (con .ge. 0.99d0) then
                     con = 0.88d0
                 endif

                if (con.lt.cut) then
                    kcut   = 1
                    biprnt = 7
                    cut = con
                endif
             endif
         elseif (A_KeyCreep.eq.2 .and. A_CrLimit.gt.0.0d0) then
c --- implicit creep
             crcrit = A_CrLimit * ( 1.d0 + A_CrTol )
             if (A_CrMax .gt. crcrit ) then
                 con = A_CrLimit / A_CrMax
c --- same principle as plasticity
                 if (con .ge. 0.99d0) then
                     con = 0.88d0
                 endif
                 if (con .lt. cut) then
                     kcut   = 1
                     biprnt = 7
                     cut = con
                 endif
            endif
         endif
c
c --- for displacement limit then do A_KeyBisect
         if (abs(A_MaxDofInc).gt.abs(A_DispLimit) .and. kcut.ne.1) then
             kcut = 1
             biprnt = 8
             cut = 0.5d0
         endif
      endif

c --- set next time increment and printout
      if (kcut.ge.1) then
         if (mintime .and. A_KeyStab.gt.0
     &       .and. A_SubStepNum.eq.1) then
             cut = 1.0d0
             kStbCut = 1
             A_StabFirstStep = 3
         else
             kStbCut = 0
         endif

         if (mintime) then
c --- reset cut for the case of first bisection and minimum time step
             cut = 1.0d0
         elseif (C_nConA.gt.0 .and. A_BisecNum.eq.0) then
c --- if some contacts were deleted, we do not decrease the time step
             cut = 1.0d0
         endif    

         con = cut*A_TimeInc
         if (con .ge. 1.5d0*A_MinTimeStep) then
c --- use half time inc
             dtnext = con
         else
c --- use min time inc
             dtnext = A_MinTimeStep
         endif
c --- bisection counter
         A_BisecNum = A_BisecNum + 1
         
c --- reset time
         A_TimeCur = A_TimeCur - A_TimeInc
c --- new time increment
         if (biloc.eq.8) A_TimeIno = A_TimeInc

         A_TimeInc = dtnext

         if (biprnt.eq.1) then
            k = A_IterMax
            if (newcutback.eq.0) then
               write (iott,3001) C_nIterEsti,k
            else
               write (iott,3091) A_IterNum
            endif
 3001       format (4x,'--> PREDICTED NUMBER OF ITERATIONS FOR ',
     &       'CONVERGENCE=',i4,' EXCEEDS MAX OF',i4)
 3091       format (4x,'--> SOLUTION PATTERN DIVERGES ',
     &       'AT ITERATION =',i4)
         elseif (biprnt .eq. 2) then
            write (iott,3002)
 3002       format (4x,'--> NEGATIVE PIVOT ENCOUNTERED')
         elseif (biprnt .eq. 3) then
            write (iott,3003) A_IterNum
 3003       format (4x,'--> SOLUTION NOT CONVERGED AFTER',
     &              i4,' ITERATIONS'/)
         elseif (biprnt .eq. 5) then
            call DistortMessage(1)
         elseif (biprnt.eq.6 .or. biprnt.eq.7 .or. biprnt.eq.8) then
         call erhandler('plcrdiag',5000,1,
     &   ' Incremental plastic/creep strain is too large',0.0d0,' ')
         elseif (biprnt .eq. 9) then
            write (iott,3029) A_IterNum
 3029       format (4x,'--> SOLUTION NOT CONVERGED AFTER',
     &              i4,' ITERATIONS DUE TO CONTACT STATUS CHANGE')
         endif
         write (iott,2030) A_SubStepNum
 2030    format (' --> SUBSTEP',i6,'  NOT COMPLETED.'/)
         if (kStbCut.eq.0) then
            write (iott,3009) A_BisecNum,dtnext
 3009       format (' --> START BISECTION NUMBER',i4,'    NEW LOAD ',
     &      'INC=',g13.5)
         elseif (kStbCut.eq.1) then
            write (iott,3010) A_BisecNum,dtnext
 3010       format (' --> BEGIN BISECTION NUMBER',i4,'        TIME ',
     &      'INCREMENT=',g13.5/' --> LOAD INC REACHED MINIMUM.'
     &      ,' STABILIZATION IS ACTIVATED.')
         endif
      endif

c --- copy OVR to SVR or initialize dof vectors
      if(kcut.eq.1 .and. A_LoadStepNum.eq.1 .and. A_SubStepNum.eq.1)then
c --- initialize (first step)
         call v_zero (M_SoluVects(1),A_NumSolVects*D_ULen)
         if (M_PtrERot.ne.PTRFTNNULL)
     &       call v_zero (M_ERot(1,1),3*D_BacLen)
         A_KeyFirstPass = 1
      elseif (kcut .eq. 1) then
         call OvrToSvr
      endif


#if defined(DEBUG)
      call SubEnd('IfBisect')
#endif

      return
      end

