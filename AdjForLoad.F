      subroutine AdjForLoad
c --- primary function: modify M_DofBits vector for forces
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "dof_info.inc"
c
      integer   forget
      integer   n, i, node, ii, j
      double precision dvalue(4)

#if defined(DEBUG)
      call SubBeg('AdjForLoad')
#endif

c --- set the forces for the compressed nodes
      n = 0
      do i = 1,D_BacLen
c --- get expanded node number
         node = M_BackList(i)
c --- get the bit pattern of specified forces
          do j = 1,6
             n = n + 1
             ii = forget(node,j,dvalue(1))
             if (ii.gt.0) then
                 M_DofBits(n) = ibset(M_DofBits(n),W_ForcBit)
                 M_DofBits(n) = ibset(M_DofBits(n),W_TranBit)
             endif
          enddo
      enddo

#if defined(DEBUG)
      call SubEnd('AdjForLoad')
#endif

      return
      end
