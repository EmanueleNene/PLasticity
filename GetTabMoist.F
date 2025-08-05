      subroutine GetTabMoist (tab,xyz,dat)
c --- primary function: get the current moisture from a 4D table
c
c  input arguments:
c     tab    (int,sc,dp)           - pointer to the table
c     xyz   (dp,ar(3),in)          - coordinates
c
c  output arguments:
c     dat    (dp,sc,out)           - data at this time point
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
#include "ParTab.inc"
c
      double precision ParTabEval
      
      double precision
     & tab,xyz(3),dat
      
      double precision
     & vect(PARTABMAX)

#if defined(DEBUG)
      call SubBeg('GetTabMoist')
#endif
      vect(1) = A_TimeCur
      call v_move(xyz(1),vect(2),3)
c      call lockSMP(W_LockTab)
      call ParTabPut (vect(1))
      dat = ParTabEval(tab)
c      call unlockSMP(W_LockTab)

#if defined(DEBUG)
      call SubEnd('GetTabMoist')
#endif

      return
      end