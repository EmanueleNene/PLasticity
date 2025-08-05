      subroutine m_norm (preoper,a,d,n)
c --- primary function: normalize a square matrix by pre-multipling 
c                       the matrix with inverse diagonal
c
c  input arguments:
c    preoper (int,sc,in)          = 1 pre-multiply the matrix with inverse diagonal
c                                <> 1 post-multiply the matrix with inverse diagonal 
c    a       (dp,ar(n,n),inout)  - matrix to be scaled
c    d       (dp,ar(n),inout)    - diagonal inverse
c
#include "impcom.inc"
c
      integer n,i,j,preoper
      double precision a(n,n),d(n), dt
      if (preoper.eq.1) then
          do i = 1,n
             dt = 0.d0
             do j=1,n
                dt = max ( abs(a(i,j)), dt )
             end do
             if (dt.eq.0.0d0) then
                 dt = 1.0d0
             endif
             d(i) = 1.d0/dt
          enddo
          do i = 1,n
             do j = 1,n
                a(i,j) = a(i,j)*d(i)
             enddo
          enddo
      else
          do i = 1,n
             do j = 1,n
                a(i,j) = a(i,j)*d(j)
             enddo
          enddo
      endif

      return
      end