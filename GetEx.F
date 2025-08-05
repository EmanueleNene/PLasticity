      subroutine GetEx (elem,ebeam,matsld,mattab)
c --- primary function: retrive attached element EX
c
c  input arguments:
c       elem   (int,sc,in)         - element number
c       matsld (int,sc,in)         - element material ID
c
c  output arguments:
c       ebeam (dp,sc,out)          - Young's modulus
c       mattab (dp,sc,out)         - material behavior code
c
#include "impcom.inc"
#include "constants.inc"
c
      integer elem,mattab
      double precision ebeam
      integer matsld
      integer*8 nlkword(W_nMatCol,W_nMatRow), keyplL, keycrL

#if defined(DEBUG)
      call SubBeg('GetEx')
#endif

c --- find Young's modulus from attached beam
      call GetMatProps (elem,matsld,1,0.0d0,ebeam,1)
      call GetMatKeys (matsld,nlkword(1,1))
c --- return inelasticity key
      keyplL = nlkword(1,1)
      keycrL = nlkword(1,2)
      mattab = 0
      if ((keycrL+keyplL).gt.0) mattab = 1

#if defined(DEBUG)
      call SubEnd('GetEx')
#endif

      return
      end
