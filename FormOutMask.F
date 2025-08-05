      subroutine FormOutMask (entmsk,leng,nelem,nnodes)
c --- primary function: set output mask for element level output
c
c  input arguments:
c     leng     (int,sc,in)             - dimensioned length of entmsk (max(nelem,nnodes)
c     nelem    (int,sc,in)             - highest element number defined
c     nnodes   (int,sc,in)             - number of nodes in compressed set
c
c  output arguments:
c     entmsk   (int,ar(leng),inout)    - updated output mask
c
#include "impcom.inc"
#include "constants.inc"
#include "RST_file.inc"
#include "solu_info.inc"
c
      integer leng,entmsk(leng),nelem,nnodes,item,j,kk,cod
      logical IfEqvDP
      
#if defined(DEBUG)
      call SubBeg('FormOutMask')
#endif

      if (R_MaxNR.eq.-1) then
          cod = 101
      else
          cod = 102
      endif

      R_OutCode = 0
c --- check if output should be done at this step
      if (abs(R_NRes).gt.0 .and. R_CurSet.lt.R_MaxSet .and.
     &    R_CurSet.gt.0 .and. R_Stop.eq.0) then
          if (R_NRes.gt.0) then
c --- save cirtain number of steps
              if (AINT(R_NRes*A_CurLoadFact)+1 .gt. R_CurSet) then
                  R_OutCode = cod
              endif
          else
c --- save every (R_NRes)'th step
              if (mod(A_SubStepNum,R_NRes).eq.0) then
                  R_OutCode = cod
              endif
          endif
      elseif (A_BisecNum.gt.0 .or. A_TimeInc.le.A_MinTimeStep) then
          R_OutCode = cod
      elseif (R_Stop.eq.-1) then
          R_OutCode = cod
      else
c --- first and last steps are always saved (R_Stop=-2 will go here)
          R_OutCode = cod
      endif
      
c --- adjust output mask for each specification
      if (R_OutCode.gt.0) then
          kk = max(nnodes,nelem)
          do item = 1,kk
             call SetOutMask (R_OutCode,item,nnodes,nelem,entmsk(item))
          enddo
      endif

c --- define global output mask
      R_OutMask = 0
      do j = 1,leng
         R_OutMask = ior(R_OutMask,entmsk(j))
      enddo

#if defined(DEBUG)
      call SubEnd('FormOutMask')
#endif

      return
      end
