      subroutine m_move (a,b,na,nb,nrow,ncol)
c --- primary function: move one matrix into another matrix
c
c  input arguments:
c     a        (dp,ar(na,*),in)  - matrix a
c     na       (int,sc,in)       - number of rows (first dimension) of matrix a
c     nb       (int,sc,in)       - number of rows (first dimension) of matrix b
c     nrow     (int,sc,in)       - number of rows to exchange
c     ncol     (int,sc,in)       - number of columns to exchange
c
c  output arguments:
c     b        (dp,ar(nb,*),out) - matrix b
c
#include "impcom.inc"
c
      integer  na, nb, nrow, ncol, i, j
      double precision   a(na,*),b(nb,*)
 

      if (na.eq.nrow .and. nb.eq.nrow) then

          call v_move(a(1,1),b(1,1),nrow*ncol)

      else

          do j = 1,ncol
             do i = 1, nrow
                b(i,j) = a(i,j)
             enddo
          enddo

      endif

      return
      end
