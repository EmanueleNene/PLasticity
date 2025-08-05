      subroutine SaveReact (gmask,f,locfpL,nrf)
c --- primary function: write out the reaction forces
c
c  input arguments:
c     gmask    (int,sc,in)       - output mask
c     f        (dp,ar(*),in)     - reaction force vector
c     locfpL   (intL,sc,in)      - starting location for reaction forces
c
c  output arguments:
c     nrf      (int,sc,out)      - number of reaction forces output
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "RST_file.inc"
#include "dof_info.inc"
c
      PTRFTN    MemAlloc

      integer*8          locfpL
      integer            gmask, k,
     &                   nrf, i, j, node, ka
      double precision   f(*)

      pointer (piWorkL,iWork)
      pointer (pdWorkL,dWork)
      integer  iWork(*)
      double precision  dWork(*)

#if defined(DEBUG)
      call SubBeg('SaveReact')
#endif

c --- write out the reaction forces
      if (btest(gmask,W_OutReact)) then
          k = D_BacLen*D_NumDof
          piWorkL = MemAlloc(k,W_LenInt,'piWorkL         ')
          pdWorkL = MemAlloc(k,W_LenDouble,'pdWorkL         ')
          nrf = 0
          do i = 1,D_FwdLen
             node = M_FrwdList(i)
             if (node .le. 0) cycle
             ka = (node-1)*D_NumDof
             do j = 1,D_NumDof
                if (btest(M_DofBits(ka+j),W_DispBit)) then ! constraint
                    nrf = nrf + 1
                    dWork(nrf) = f(ka+j)
                    iWork(nrf) = (node-1)*D_NumDof + j
                endif
             enddo
          enddo
          if (nrf .gt. 0) then
              call Buf_Write (R_BlockNum,locfpL,nrf,iWork(1),nrf)
              call Buf_Write (R_BlockNum,locfpL,nrf*W_IntPerDp,
     &                        dWork(1),0)
          endif
          call MemFree (pdWorkL)
          call MemFree (piWorkL)
      endif


#if defined(DEBUG)
      call SubEnd('SaveReact')
#endif

      return
      end
