      subroutine UpdTransParams
c --- primary function: calculate transient integration constants
c
#include "impcom.inc"
#include "solu_info.inc"
#include "dyn_info.inc"
c
      double precision  alpha,delta,con,alpha3
      equivalence (alpha,I_Consts(2)),(delta,I_Consts(3))

c --- only implicit now
      I_Opt = 0
      I_Consts(1) = 0.005d0
      I_Consts(2) = 0.25d0*(1.0d0 + I_Consts(1))**2
      I_Consts(3) = 0.5d0 + I_Consts(1)

      I_RespStiff = 0.0d0
      I_RespMass = 0.0d0
      if (I_Opt.eq.0) then
c --- 2nd order implicit, Newmark
         I_Prm(9) = 0.0d0
         I_Prm(10)= 0.0d0
         I_Prm(11)= 1.0d0
         I_Prm(12)= 1.0d0
         I_Prm(1) = 1.0d0/(alpha*A_TimeInc*A_TimeInc)
         I_Prm(2) = delta/(alpha*A_TimeInc)
         I_Prm(3) = (0.5d0-alpha)/alpha
         I_Prm(4) = (0.5d0*delta-alpha)*A_TimeInc/alpha
         I_Prm(5) = 1.0d0/(alpha*A_TimeInc)
         I_Prm(6) = (delta-alpha)/alpha
         I_Prm(7) = (1.0d0-delta)*A_TimeInc
         I_Prm(8) = delta*A_TimeInc
         alpha3 = delta - alpha - 0.5d0
         I_Prm(13)= 1.0d0
         I_Prm(14)= 0.0d0
         I_Prm(15)= 0.0d0
         I_Prm(16)= 1.0d0
         I_Prm(17)= 0.0d0
         I_Prm(18)= 0.0d0
         I_Prm(19)= 1.0d0/(alpha*A_TimeInc*A_TimeInc)
         I_Prm(20)= 1.0d0/(alpha*A_TimeInc)
         I_Prm(21)= 0.5d0/alpha - 1.0d0
      elseif (I_Opt.eq.2) then
c --- 2nd order explicit
         call v_zero (I_Prm(1),8)
         con = A_TimeIno + A_TimeInc
c --- protect against zero divide
         if (con .gt. 0.0d0) then
             if (A_TimeInc .gt. 0.0d0) then
                 I_Prm(1) = 2.0d0/(A_TimeInc*con)
                 I_Prm(2) = A_TimeIno/(A_TimeInc*con)
             endif
             if (A_TimeIno .gt. 0.0d0) then
                 I_Prm(3) = 2.0d0/(A_TimeIno*con)
                 I_Prm(4) = -A_TimeInc/(A_TimeIno*con)
             endif
         endif
      endif
      alpha3 = delta - alpha - 0.5d0
      I_Prm(22) = 1.0d0/alpha3/A_TimeInc/A_TimeInc
      I_Prm(23) =-1.0d0/alpha3/A_TimeInc
      I_Prm(24) = (delta-alpha)/alpha3
      I_Prm(25) =-(1.0d0-delta)/alpha3/A_TimeInc
      I_Prm(26) = (0.5d0-alpha)/alpha3
      I_Prm(27) =-(0.5d0*delta-alpha)/alpha3*A_TimeInc
      I_Prm(28) = 1.0d0 - delta
      I_Prm(29) = (0.5d0-alpha)*A_TimeInc*A_TimeInc
      I_Prm(30) = alpha*A_TimeInc*A_TimeInc
      I_Prm(31) = delta/alpha/A_TimeInc
      I_Prm(32) = (0.5d0*delta/alpha-1.0d0)*A_TimeInc
      I_Prm(33) = delta/alpha-1.0d0

      return
      end
