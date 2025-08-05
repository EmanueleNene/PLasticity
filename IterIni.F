      subroutine IterIni (dincmx)
c --- primary function: finish up an equilibium iteration
c
c  input arguments:
c     dincmx   (dp,sc,in)          - maximum displacement increment
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
c
      integer  GetWriteUnit
      logical  success,commitqq
      integer  iott
      double precision  dincmx

#if defined(DEBUG)
      call SubBeg('IterIni')
#endif

      iott = GetWriteUnit()

c --- write out iteration trailer
      if (A_IterMax.gt.1) then
          write (iott,2010) A_IterNum,dincmx
 2010     format (4x,'ITER',i4,' FINISHED.  ',
     &          ' MAX DISPLACEMENT INC = ',g11.4)
      endif

c --- flush the buffer
#if defined (PCWINNT_SYS) || defined(LINUX_SYS)
      success = commitqq (iott)
#else
      call flush(iott)
#endif

#if defined(DEBUG)
      call SubEnd('IterIni')
#endif

      return
      end
