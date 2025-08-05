      subroutine AdjForBCs
c --- primary function: modify M_DofBits(*) vector for displacements
c
#include "impcom.inc"
#include "constants.inc"
#include "dof_info.inc"
#include "in_mem.inc"
c
      integer  conchk, disget

      integer  n,i,node,j,ii
      double precision dval(4)
      
#if defined(DEBUG)
      call SubBeg('AdjForBCs')
#endif

c --- set the displacements for the compressed nodes
      n = 0
      do i = 1,D_BacLen
c --- get expanded node number
         node = M_BackList(i)
c --- get the bit pattern of specified displacements
         do j = 1,6
            n = n + 1
            if (conchk(node,j).ne.2) then
                M_DofBits(n) = ibset(M_DofBits(n),W_DispBit)
                ii = disget(node,j,dval(1))
                if (i.gt.0)
     &              M_DofBits(n) = ibset(M_DofBits(n),W_TranBit)
            endif
         enddo
      enddo

#if defined(DEBUG)
      call SubEnd('AdjForBCs')
#endif

      return
      end


