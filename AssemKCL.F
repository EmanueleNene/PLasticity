      subroutine AssemKCL (plSpL,nr,kelfil,ls,Kmtrx)
c --- primary function: assemble elements
c
c  input arguments:
c     plSpL    (PTR,sc,in)       - pointer to sparse matrix
c     nr       (int,sc,in)       - number of rows in element matrix
c     kelfil   (int,ar(10),in)   - flags for existance of matrices
c     ls       (int,ar(nr),in)   - the degrees of freedom
c     Kmtrx    (dp,ar(nr,nr),in) - the element matrix
c
#include "impcom.inc"
#include "dof_info.inc"
c

      PTRFTN             plSpL
      integer            nr, kelfil(10), ls(nr)
      double precision   Kmtrx(nr,nr,*), dummy


#if defined(DEBUG)
      call SubBeg('AssemKCL')
#endif

      if (kelfil(1) .gt. 0) then
          if (D_NumCp .gt. 0) call ConvertDofCp (nr,ls(1))
c --- add to sparse matrix
          call SpAddMatrix (plSpL,nr,Kmtrx(1,1,1),ls(1))
      endif

#if defined(DEBUG)
      call SubEnd('AssemKCL')
#endif

      return
      end
