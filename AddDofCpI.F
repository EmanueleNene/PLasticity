      subroutine AddDofCpI (iLeng,iVec)
c --- primary function: add slave DOF to input vector
c
c  input arguments:
c     iLeng    (int,sc,in)         - length of vector to be expanded
c
c  output arguments:
c     iVec     (int,ar(*),inout)   - expanded vector
c
#include "impcom.inc"
#include "dof_info.inc"
#include "in_mem.inc"
c
      integer   i, k, iLeng, leng, iVec(*)


#if defined(DEBUG)
      call SubBeg('AddDofCpI')
#endif

c --- get the coupled DOF vector
      if (D_NumCp .eq. 0) goto 999
      leng = iLeng
      if (D_ULen .lt. iLeng) leng = D_ULen

c --- loop through the coupled DOF vector
      do i = 1,leng
         k = M_SlaveCpDofs(i)
         if (k .gt. 0) iVec(i) = iVec(k)
      enddo

  999 continue

#if defined(DEBUG)
      call SubEnd('AddDofCpI')
#endif

      return
      end
