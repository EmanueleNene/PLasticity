      subroutine ConvertDofCp (iLeng,iVec)
c --- primary function: convert DOF vector to master DOFs
c
c  input arguments:
c     iLeng    (int,sc,in)          - length of vector to be converted
c
c  output arguments:
c     iVec     (int,ar(*),inout)    - converted vector
c
#include "impcom.inc"
#include "dof_info.inc"
#include "in_mem.inc"
c
      integer  iLeng, iVec(*), i, leng, k

c --- get the coupled DOF vector
      if (D_NumCp .eq. 0) goto 999
      leng = iLeng
      if (D_ULen .lt. iLeng) leng = D_ULen

c --- loop through the coupled DOF vector
      do i = 1,leng
         if (iVec(i) .le. 0) cycle
         k = M_SlaveCpDofs(iVec(i))
         if (k .gt. 0) iVec(i) = k
      enddo

  999 continue

      return
      end
