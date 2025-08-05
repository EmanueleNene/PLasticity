      function m_solve (a,nd,n,nc)
c --- primary function:     solve a set of simultaneous linear equations 
c
c --- secondary functions:  invert a matrix (any)
c
c  input arguments:
c     a        (dp,ar(nd,*),inout)  - matrix to be solved or inverted
c                                     second dimension must be at least: n + abs(nc)
c     nd       (int,sc,in)          - first dimension of the a matrix
c     n        (int,sc,in)          - number of equations
c     nc       (int,sc,in)          - number of additional columns.  
c                                     * if nc = +n or -n, invert n x n matrix and 
c                                     put result in the n+1 to 2xn columns.
c                                     * if nc is 0 or negative, nc will be reset to 
c                                     n and then m_solve will set up identity 
c                                     matrix after the input matrix, where the 
c                                     result of the inversion will be put.
c                                     * if nc is positive and less than n, do a
c                                     partial inversion. see example 1 below.
c                                             
c  output arguments:
c     m_solve  (in,sc,out)         - 0 - non-singular matrix
c                                    1 - singular matrix
c                                    2 - near-singular matrix
c     a        (dp,ar(nd,*),inout) - results or inverted matrix.  
c                                    starts in column n+1.
c                                    note:  original information is destroyed.
c
c  example 1:   solve three simultaneous linear equations:
c               i = m_solve (a(1,1),3,3,1)
c               calling routine has a dimensioned as a(3,4)
c               each equation has its 3 coefficents in the first 3 columns,
c               and the constant term is in the fourth column.
c               solution is in fourth column.
c
c  example 2:   invert a 3x3 matrix:
c               i = m_solve (a(1,1),3,3,-3)
c               calling routine has a dimensioned as a(3,6)
c               input matrix was input in first 3 columns
c               output matrix in ouput in last 3 columns
c           
#include "impcom.inc"
c                     
      integer  m_solve, nd, n, nc
      double precision  a(nd,*)
                                      
      integer  n1, npn1, i, j, k, ip1
      double precision  con, TINYVAL, tol, minDiag
      parameter        (TINYVAL=1.0d-16)

      m_solve = 0
      if (n .le. 0) goto 999
      n1 = nc
             
c --- set up for matrix inverse
      if (nc .le. 0) then
          n1 = n
          k = n
          do i = 1,n
             k = k + 1
             do j = 1,n
                a(j,k) = 0.0d0
             enddo
          enddo 
          k = n
          do i = 1,n
             k = k + 1
             a(i,k) = 1.0d0
          enddo
      endif
      npn1 = n + n1

c --- tol = TINY was too tight...
c --- get min and max absolute diagonal values to use to scale tolerance
      if (n .gt. 1) then
          minDiag = HUGE
          do i = 1,n
             con = abs(a(i,i))
             if (con .lt. minDiag) minDiag = con
          enddo
          if (minDiag .le. 0.0d0) minDiag = 1.0d0
          tol = minDiag * TINYVAL
      else
c --- 1x1 case is trivial, with no loss of accuracy, so let it go through
          tol = 0.0d0
      endif

c --- forward pass
      do i = 1,n
c --- test for a near-singular or singular matrix
         if (abs(a(i,i)) .le. tol) then
             m_solve = 2
             if (abs(a(i,i)) .le. 0.0d0) m_solve = 1
             goto 999
         endif

         ip1 = i + 1
         con = 1.0d0/a(i,i)
         do k = ip1,npn1
            a(i,k) = con*a(i,k)
         enddo
         do k = ip1,n
            con = -a(k,i)
!dir$ ivdep
            do j = i+1,npn1
               a(k,j) = a(k,j) + con*a(i,j)
            enddo
         enddo 
      enddo
              
c --- back substitution
      do j = n+1,npn1
         do k = n-1,1,-1
            do i = k+1,n
               a(k,j) = a(k,j) - a(i,j)*a(k,i)
            enddo
         enddo
      enddo

 999  return
      end
