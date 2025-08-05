      subroutine WriteToRST (edtype,lcerstL,edindxL,nval,dvalue)
c --- primary function: output element data to result file
c
c  input arguments:
c        edtype    (int,sc,in)        - element data type
c        lcerstL   (intL,sc,inout)    - pointer to results file position
c        edindxL   (intL,ar(25),inout)- index to results file data
c        nval      (int,sc,in)        - the total number of values
c        dvalue    (dp,ar(nval),in)   - output values
c
#include "impcom.inc"
#include "constants.inc"
#include "RST_file.inc"
c
      integer   InqSMP

      integer*8         lcerstL, edindxL(25)
      integer           edtype, nval, dtype, i, n, nprocs
      double precision  dvalue(nval)

#if defined(DEBUG)
      call SubBeg('WriteToRST')
#endif

      if (nval .le. 0) goto 999

      dtype = abs(edtype)
      n = nval

c --- check for all zero record
      nprocs = InqSMP(0)
      if (nprocs .le. 1) then
          do 10 i = 1,nval
             if (dvalue(i) .ne. 0.0d0) goto 50
  10      continue
          edindxL(dtype) = -n
          goto 999
      endif

c --- output the record
  50  edindxL(dtype) = lcerstL
      call lockSMP (W_LockRst)
      call Buf_Write (R_BlockNum,lcerstL,nval*W_IntPerDp,dvalue(1),0)
      call unlockSMP (W_LockRst)


  999 continue

#if defined(DEBUG)
      call SubEnd('WriteToRST')
#endif

      return
      end
