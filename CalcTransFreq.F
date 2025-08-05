      subroutine CalcTransFreq
c --- primary function: calculate the response frequency 
c                                                                     
#include "impcom.inc"
#include "solu_info.inc"
#include "constants.inc"
#include "dyn_info.inc"
c
      integer  i
      double precision  con
                                    
      I_RespFrq = 0.0d0
      if (A_Trans.gt.0 .and. A_TimeInc.gt.0.0d0) then
c --- second order
          if (I_RespStiff.gt.0.0d0 .and. I_RespMass.gt.0.0d0) then
c --- check for possible overflow and underflow conditions
             con = log10(I_RespStiff) - log10(I_RespMass)
             i = nint(con)
             if (i.gt.-30 .and. i.lt.30) then
                 I_RespFrq = I_RespStiff/I_RespMass
             endif
             if (I_RespFrq.gt.0.0d0) then
                 I_RespFrq = sqrt(I_RespFrq)/(2.0d0*W_Pi)
             endif
          endif
      endif

      return
      end
      

