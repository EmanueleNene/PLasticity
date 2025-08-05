      subroutine GetD (node,idf,dval)
c --- primary function: get displacements
c
c  input arguments:
c      inode    (int,sc,in)    - node number 
c      idf      (int,sc,in)    - reference number for the DOF
c
c  output arguments:
c      value    (dp,ar(4),out) - constraint values
c
#include "impcom.inc"
#include "solu_info.inc"
#if defined(MOIST)
c#include "constants.inc"
#endif
c
      integer disget
      double precision  ParTabEval
            
      integer node, idf, i
      double precision  dval(4)
      logical ParTabChek         
c --- dval should be zeroed in disget      
      i = disget (node,idf,dval(1)) 
         
      if (i.gt.0) then
          if (ParTabChek(dval(1))) then
#if defined(MOIST)
c              call lockSMP(W_LockTab)
#endif              
              call ParTabPut1 (1,A_TimeCur)
              dval(1) = ParTabEval(dval(1))
#if defined(MOIST)
c              call unlockSMP(W_LockTab)
#endif              
              if (dval(1).eq.HUGE) then
                  dval(1) = 0.0d0
              endif
          else                           
c --- compute load vector contribution
              dval(1) = dval(1)*A_CurLoadFact
          endif
      endif

      return
      end