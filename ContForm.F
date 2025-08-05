      subroutine ContForm (kerr)
c --- primary function: set up contact data base
c
c
#include "impcom.inc"
#include "constants.inc"
#include "dof_info.inc"
#include "cont_info.inc"
#include "in_mem.inc"

      integer   ndgall

      integer            inode,kloc,nodei, kerr, jrest, j
      double precision   coor(6)

#if defined(DEBUG)
      call SubBeg('ContForm')
#endif

      kerr = 0
      do inode = 1,D_BacLen
         kloc = D_NumDof*(inode-1) + 1
         nodei = M_BackList(inode)
         j = ndgall (nodei,coor(1))
         if (j .ne. 0) then
             jrest = D_NumDof - 3
             if (jrest .gt .0) then
                 call v_move (coor(1),K_NdCoord(kloc),3)
                 call v_zero (K_NdCoord(kloc+3),jrest)
             else
                 call v_move (coor(1),K_NdCoord(kloc),D_NumDof)
             endif
         else
             call v_zero (K_NdCoord(kloc),D_NumDof)
         endif
      enddo
      call ContManager (D_NumEl,D_BacLen,D_NumDof,
     &                  M_FrwdList(1),M_BackList(1),K_nPerSeg,K_nSeg,
     &                  K_nPerFounI,K_nMast,K_nPerMast,K_nPerFounDp,
     &                  K_nFoun,K_NdCoord(1),M_SoluVects(1),
     &                  K_Seg(1),K_FounI(1),K_FounDp(1),
     &                  K_ContInfo(1),K_PiloNd(1),kerr)
c

#if defined(DEBUG)
      call SubEnd('ContForm')
#endif

      return
      end



