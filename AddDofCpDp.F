      subroutine AddDofCpDp (uLeng,uVec)
c --- primary function: add slave DOF to input vector
c
c  input arguments:
c     uLeng    (int,sc,in)         - length of vector to be expanded
c
c  output arguments:
c     uVec     (dp,ar(*),inout)    - expanded vector
c
#include "impcom.inc"
#include "dof_info.inc"
#include "in_mem.inc"
c
      integer   i, k, uLeng, leng
      double precision  uVec(*)

#if defined(DEBUG)
      call SubBeg('AddDofCpDp')
#endif

c --- get the coupled DOF vector
      if (D_NumCp .eq. 0) goto 999
      leng = uLeng
      if (D_ULen .lt. uLeng) leng = D_ULen

c --- loop through the coupled DOF vector
      do i = 1,leng
         k = M_SlaveCpDofs(i)
         if (k.gt.0 .and. uVec(i).ne.HUGE) uVec(i) = uVec(k)
      enddo

  999 continue

#if defined(DEBUG)
      call SubEnd('AddDofCpDp')
#endif

      return
      end
