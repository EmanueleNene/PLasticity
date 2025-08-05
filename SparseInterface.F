      subroutine SparseInterface (plSpL)
c --- primary function: interface to sparse solver assembly
c
c  input arguments:
c    nume      (int,sc,in)       - number of elements (compressed)
c
c  output arguments:
c    plSpL     (PTR,sc,out)      - pointer to sparse matrix
c
#include "impcom.inc"
#include "constants.inc"
#include "dof_info.inc"
c
      PTRFTN    plSpL


#if defined(DEBUG)
      call SubBeg('SparseInterface')
#endif

      call SpSymbAssem (D_NumDof,plSpL)
      call SpAlloc     (plSpL,1)

#if defined(DEBUG)
      call SubEnd('SparseInterface')
#endif

      return
      end
