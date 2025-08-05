      subroutine cAtxA_C(a,c, nc,n, con)
c --- primary function: do C = con*(At*A) 
c
c  input arguments:
c     a        (dp,ar(*),in)       - vector to be multiplied by itself
c     nc       (int,sc,in)         - number of rows in the c matrix
c     n        (int,sc,in)         - size of square matrix
c     con      (dp,sc,in)          - multiplier on above square matrix
c                                             
c  output arguments:
c     c        (dp,ar(nc,*),inout) - matrix to be accumulated onto
c                                    only the lower triangular matrix is done
c
c
#include "impcom.inc"
c
      integer            nc, n, i, j, irem
      double precision   a(*), c(nc,*), con, w(1000)

      if (con .eq. 1.0d0) then
          do i = 1,n
             w(i) = a(i)
          enddo
      else
          do i = 1,n
             w(i) = a(i)*con
          enddo
      endif

c --- do the shortest vector lengths here
      if (n .gt. 8) goto 9
      goto ( 1,2,3,4,5,6,7,8 ), n

 1    c(1,1) = c(1,1) + w(1)*a(1)
      goto 999 

 2    c(1,1) = c(1,1) + w(1)*a(1)
      c(2,1) = c(2,1) + w(2)*a(1)
      c(2,2) = c(2,2) + w(2)*a(2)
      goto 999 

!dir$ ivdep,novector
 3    do i = 1,3
         c(i,1) = c(i,1) + w(i)*a(1)
         c(i,2) = c(i,2) + w(i)*a(2)
         c(i,3) = c(i,3) + w(i)*a(3)
      enddo
      goto 999 

!dir$ ivdep,novector
 4    do i = 1,4
         c(i,1) = c(i,1) + w(i)*a(1)
         c(i,2) = c(i,2) + w(i)*a(2)
         c(i,3) = c(i,3) + w(i)*a(3)
         c(i,4) = c(i,4) + w(i)*a(4)
      enddo
      goto 999 

!dir$ ivdep,novector
 5    do i = 1,5
         c(i,1) = c(i,1) + w(i)*a(1)
         c(i,2) = c(i,2) + w(i)*a(2)
         c(i,3) = c(i,3) + w(i)*a(3)
         c(i,4) = c(i,4) + w(i)*a(4)
         c(i,5) = c(i,5) + w(i)*a(5)
      enddo
      goto 999 

!dir$ ivdep,novector
 6    do i = 1,6
         c(i,1) = c(i,1) + w(i)*a(1)
         c(i,2) = c(i,2) + w(i)*a(2)
         c(i,3) = c(i,3) + w(i)*a(3)
         c(i,4) = c(i,4) + w(i)*a(4)
         c(i,5) = c(i,5) + w(i)*a(5)
         c(i,6) = c(i,6) + w(i)*a(6)
      enddo
      goto 999 

!dir$ ivdep,novector
 7    do i = 1,7
         c(i,1) = c(i,1) + w(i)*a(1)
         c(i,2) = c(i,2) + w(i)*a(2)
         c(i,3) = c(i,3) + w(i)*a(3)
         c(i,4) = c(i,4) + w(i)*a(4)
         c(i,5) = c(i,5) + w(i)*a(5)
         c(i,6) = c(i,6) + w(i)*a(6)
         c(i,7) = c(i,7) + w(i)*a(7)
      enddo
      goto 999 

!dir$ ivdep,novector
 8    do i = 1,8
         c(i,1) = c(i,1) + w(i)*a(1)
         c(i,2) = c(i,2) + w(i)*a(2)
         c(i,3) = c(i,3) + w(i)*a(3)
         c(i,4) = c(i,4) + w(i)*a(4)
         c(i,5) = c(i,5) + w(i)*a(5)
         c(i,6) = c(i,6) + w(i)*a(6)
         c(i,7) = c(i,7) + w(i)*a(7)
         c(i,8) = c(i,8) + w(i)*a(8)
      enddo
      goto 999 

c --- lengths from 9 to 64 are handled here
 9    if (n .gt. 64) goto 100
         irem = mod(n,8)
         goto (99,10,20,30,40,50,60,70), irem + 1
        
cdir$ ivdep,novector
 10      do i = 1,n
            c(i,1) = c(i,1) + w(i)*a(1)
         enddo
         goto 99

cdir$ ivdep,novector
 20      do i = 1,n
            c(i,1) = c(i,1) + w(i)*a(1)
            c(i,2) = c(i,2) + w(i)*a(2)
         enddo
         goto 99

cdir$ ivdep,novector
 30      do i = 1,n
            c(i,1) = c(i,1) + w(i)*a(1)
            c(i,2) = c(i,2) + w(i)*a(2)
            c(i,3) = c(i,3) + w(i)*a(3)
         enddo
         goto 99

cdir$ ivdep,novector
 40      do i = 1,n
            c(i,1) = c(i,1) + w(i)*a(1)
            c(i,2) = c(i,2) + w(i)*a(2)
            c(i,3) = c(i,3) + w(i)*a(3)
            c(i,4) = c(i,4) + w(i)*a(4)
         enddo
         goto 99

cdir$ ivdep,novector
 50     do i = 1,n
            c(i,1) = c(i,1) + w(i)*a(1)
            c(i,2) = c(i,2) + w(i)*a(2)
            c(i,3) = c(i,3) + w(i)*a(3)
            c(i,4) = c(i,4) + w(i)*a(4)
            c(i,5) = c(i,5) + w(i)*a(5)
         enddo
         goto 99

cdir$ ivdep,novector
 60      do i = 1,n
            c(i,1) = c(i,1) + w(i)*a(1)
            c(i,2) = c(i,2) + w(i)*a(2)
            c(i,3) = c(i,3) + w(i)*a(3)
            c(i,4) = c(i,4) + w(i)*a(4)
            c(i,5) = c(i,5) + w(i)*a(5)
            c(i,6) = c(i,6) + w(i)*a(6)
         enddo
         goto 99

cdir$ ivdep,novector
 70      do i = 1,n
            c(i,1) = c(i,1) + w(i)*a(1)
            c(i,2) = c(i,2) + w(i)*a(2)
            c(i,3) = c(i,3) + w(i)*a(3)
            c(i,4) = c(i,4) + w(i)*a(4)
            c(i,5) = c(i,5) + w(i)*a(5)
            c(i,6) = c(i,6) + w(i)*a(6)
            c(i,7) = c(i,7) + w(i)*a(7)
         enddo

 99   do j = irem+1,n,8
cdir$ ivdep,novector
         do i = j,n
            c(i,j  ) = c(i,j  ) + w(i)*a(j  )
            c(i,j+1) = c(i,j+1) + w(i)*a(j+1)
            c(i,j+2) = c(i,j+2) + w(i)*a(j+2)
            c(i,j+3) = c(i,j+3) + w(i)*a(j+3)
            c(i,j+4) = c(i,j+4) + w(i)*a(j+4)
            c(i,j+5) = c(i,j+5) + w(i)*a(j+5)
            c(i,j+6) = c(i,j+6) + w(i)*a(j+6)
            c(i,j+7) = c(i,j+7) + w(i)*a(j+7)
         enddo
      enddo
      goto 999 

c --- lengths > 64 handled here (leave out the novector directive)
 100  irem = mod(n,8)
      goto (199,110,120,130,140,150,160,170), irem + 1
        
cdir$ ivdep
 110  do i = 1,n
         c(i,1) = c(i,1) + w(i)*a(1)
      enddo
      goto 199

cdir$ ivdep
 120  do i = 1,n
         c(i,1) = c(i,1) + w(i)*a(1)
         c(i,2) = c(i,2) + w(i)*a(2)
      enddo
      goto 199

cdir$ ivdep
 130  do i = 1,n
         c(i,1) = c(i,1) + w(i)*a(1)
         c(i,2) = c(i,2) + w(i)*a(2)
         c(i,3) = c(i,3) + w(i)*a(3)
      enddo
      goto 199

cdir$ ivdep
 140  do i = 1,n
         c(i,1) = c(i,1) + w(i)*a(1)
         c(i,2) = c(i,2) + w(i)*a(2)
         c(i,3) = c(i,3) + w(i)*a(3)
         c(i,4) = c(i,4) + w(i)*a(4)
      enddo
      goto 199

cdir$ ivdep
 150  do i = 1,n
         c(i,1) = c(i,1) + w(i)*a(1)
         c(i,2) = c(i,2) + w(i)*a(2)
         c(i,3) = c(i,3) + w(i)*a(3)
         c(i,4) = c(i,4) + w(i)*a(4)
         c(i,5) = c(i,5) + w(i)*a(5)
      enddo
      goto 199

cdir$ ivdep
 160  do i = 1,n
         c(i,1) = c(i,1) + w(i)*a(1)
         c(i,2) = c(i,2) + w(i)*a(2)
         c(i,3) = c(i,3) + w(i)*a(3)
         c(i,4) = c(i,4) + w(i)*a(4)
         c(i,5) = c(i,5) + w(i)*a(5)
         c(i,6) = c(i,6) + w(i)*a(6)
      enddo
      goto 199

cdir$ ivdep
 170  do i = 1,n
         c(i,1) = c(i,1) + w(i)*a(1)
         c(i,2) = c(i,2) + w(i)*a(2)
         c(i,3) = c(i,3) + w(i)*a(3)
         c(i,4) = c(i,4) + w(i)*a(4)
         c(i,5) = c(i,5) + w(i)*a(5)
         c(i,6) = c(i,6) + w(i)*a(6)
         c(i,7) = c(i,7) + w(i)*a(7)
      enddo

 199  do j = irem+1,n,8
cdir$ ivdep
         do i = j,n
            c(i,j  ) = c(i,j  ) + w(i)*a(j  )
            c(i,j+1) = c(i,j+1) + w(i)*a(j+1)
            c(i,j+2) = c(i,j+2) + w(i)*a(j+2)
            c(i,j+3) = c(i,j+3) + w(i)*a(j+3)
            c(i,j+4) = c(i,j+4) + w(i)*a(j+4)
            c(i,j+5) = c(i,j+5) + w(i)*a(j+5)
            c(i,j+6) = c(i,j+6) + w(i)*a(j+6)
            c(i,j+7) = c(i,j+7) + w(i)*a(j+7)
         enddo
      enddo

 999  return 
      end
