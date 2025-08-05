      subroutine SumForce (force)
c --- primary function: add applied nodal forces to total force vector
c
c  output arguments:
c     force    (dp,ar(*),inout)  - global force vector
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "dof_info.inc"
#include "solu_info.inc"
c
      integer   foriqr
      integer   i,node,idf, kerr
      double precision  force(*), dvalue(4)

#if defined(DEBUG)
      call SubBeg('SumForce')
#endif

c --- add forces to load vector
      if (foriqr(0,DB_MAXDEFINED).gt.0) then
          do i = 1,D_ULen
c --- skip slave node on coupled dof (forces moved to master)
             if (btest(M_DofBits(i),W_CpBit)) cycle
             if (btest(M_DofBits(i),W_ForcBit)) then
                 call InqDof (i,node,idf)
                 call GetF (node,idf,dvalue(1))
                 force(i) = force(i) + dvalue(1)
             endif
          enddo
      endif

#if defined(DEBUG)
      call SubEnd('SumForce')
#endif

      return
      end
