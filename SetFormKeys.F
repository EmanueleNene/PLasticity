      subroutine SetFormKeys
c --- primary function: set the request keys for element matrix formation
c
c  output arguments:
c     A_MatrReq (int,ar(10),out)      - matrix and load vector form requests
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
#include "conv_stat.inc"
c
      integer kloc(4), cnvfld(4), i, k, j, j1
      data    kloc / 1, 3, 2, 4/
      double precision factor, cof, cofd
      logical IfEqvDP, mintime

#if defined(DEBUG)
      call SubBeg('SetFormKeys')
#endif
      
      mintime = ifEqvDP(A_TimeInc,A_MinTimeStep) .or. 
     &                 (A_TimeInc.lt.A_MinTimeStep)
      
      
c --- determine global options from the currently active DOF options (see ConIni.F)
      cnvfld (1) = 1         ! F
      cnvfld (3) = 1         ! M
      if (A_SpecVer.gt.0) then
          cnvfld (2) = 2     ! U
          if (A_SpecVer.eq.3 .and. C_nPivots.eq.0) cnvfld(2) = 0     ! U
          if (A_SpecVer.eq.4) cnvfld(2) = 0                          ! U
          cnvfld (4) = 0     ! R
      else
          cnvfld (2) = 0     ! U
          cnvfld (4) = 0     ! R
      endif

c --- reset C_KeyMono to ignore countings before C_IterCtrl0, the loweset is 5
      if (A_IterNum.eq.C_IterCtrl0) C_KeyMono = 0

c --- initialize the element keys
      call KeysIni
c --- set the element formation request keys
c --- k matrix
      A_MatrReq(1) = 1
c --- m matrix
      A_MatrReq(2) = 1
c --- c matrix
      A_MatrReq(3) = 1
c --- ss matrix
      A_MatrReq(4) = 1

c --- applied f vector
      A_MatrReq(5) = 1
c --- restoring force f vector
      A_MatrReq(6) = 1

c --- set newton-raphson updating keys
c --- never update history variables - only in the output pass
      A_KeyMatUpd = 0

c --- element boundary condition updating keys - do not update in form pass
      A_KeyBCsUpd = 0
      A_KeyMoistUpd = 0

c --- first iteration key
      if (A_IterNum.eq.1) A_KeyFirstIter = 1

c --- loosen up convergence criterion (this will save a few valuable iterations)
      cofd = 0.88d0 
      if (cnvfld(2).gt.0) then
          cofd = 0.6d0
      else
          if (mintime .or. A_BisecNum.gt.1) then
              cofd = 0.6d0
          endif
      endif

      do k = 1,4
         i = kloc(k)
         if (A_BisecNum.ge.1 .or. mintime) then
             if (k.le.2) then
                 cof = dble(A_BisecNum+1)
             else
                 cof = 1.1d0**A_BisecNum
             endif
         else
             cof = 1.0d0
         endif
         j1 = abs(A_IterNum)
         if (j1.eq.0) j1 = 1
         if (j1.gt.5) j1 = 5
         factor = cof/(cofd**j1)
         A_ConvTol(i) = 0.005d0*factor*dble(cnvfld(i))
      enddo
      
#if defined(DEBUG)
      call SubEnd('SetFormKeys')
#endif

      return
      end
