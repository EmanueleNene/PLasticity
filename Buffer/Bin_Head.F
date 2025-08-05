      subroutine Bin_Head (nblk,nunit,filposL)
c --- primary function:    put standard header on a RST file
c --- secondary functions: return the first data position
c
c  input arguments:
c     nblk     (int,sc,in)       - block number of open binary file
c     nunit    (int,sc,in)       - the unit number for this file
c
c  output arguments:
c     filposL  (intL,sc,out)     - the position after the header
c
c
#include "impcom.inc"
#include "constants.inc"
#include "file_sys.inc"
c
      integer  InqSMP

      integer*8  filposL
      integer  nblk, nunit, ivect(100)
      double precision  rev
      character*4  revlab


#if defined(DEBUG)
      call SubBeg('Bin_Head')
#endif

      if (nblk .le. 0) goto 999

      call vi_zero (ivect(1),100)

c --- unit number and internal/external key
      ivect(1) = nunit
      ivect(2) = -1

c --- version identifiers
      rev = 11.0d0
      write (revlab,'(F4.1)') rev
      call StringEncode (revlab,ivect(10),1)

c --- number of processors
      ivect(29) = InqSMP(0)

c --- split point of file (very big, otherwise = reading problem)
#if defined(FLIMIT)
      ivect(95) = 65536
#endif
c --- write at start of file
      filposL = 0
      B_CurAddr(nblk) = 0
      call Buf_Write (nblk,filposL,100,ivect(1),100)

  999 continue

#if defined(DEBUG)
      call SubEnd('Bin_Head')
#endif

      return
      end
