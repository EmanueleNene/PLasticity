      subroutine m_symLtoU (a,nd,n)
c --- primary function: fill upper triangle from lower triangle  
c
c  input arguments:
c     nd       (int,sc,in)         - number of rows of the a matrix
c     n        (int,sc,in)         - size of matrix to be processed
c
c  output arguments:
c     a        (dp,ar(nd,*),inout) - matrix that has its lower triangular part
c                                    copied to its upper triangular part
c
#include "impcom.inc"
c
      integer  nd, n, i, j
      double precision  a(nd,*)

      do j = 1,n-1
c --- neglect vector dependencies
!dir$ ivdep
         do i = j+1,n
            a(j,i) = a(i,j)
         enddo
      enddo

      return
      end
