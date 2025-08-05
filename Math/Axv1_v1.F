      subroutine Axv1_v1 (a,v,nr,nc)
c --- primary function: multiply vector v by matrix a
c
c  input arguments:
c     a        (dp,ar(nr,nc),in) - matrix a
c     nr       (int,sc,in)       - number of rows in matrix a
c     nc       (int,sc,in)       - number of columns to multiply in matrix a
c
c  output arguments:
c     v        (dp,ar(nr),inout) - product, stored in vector v
c
#include "impcom.inc"
c
      integer  nr, nc
      double precision  a(nr,*), v(nr), w(60)
      integer  i, j, k, n, irem

      n = abs(nc)

      irem = mod(n,8)
      goto (80,10,20,30,40,50,60,70) irem+1

!dir$ novector
 10   do i=1,nr
         w(i) = 
     +           a(i,1) * v(1)
      enddo
      goto 99

!dir$ novector
 20   do i=1,nr
         w(i) = 
     +           a(i,1) * v(1)
     +         + a(i,2) * v(2)
      enddo
      goto 99

!dir$ novector
 30   do i=1,nr
         w(i) = 
     +           a(i,1) * v(1)
     +         + a(i,2) * v(2)
     +         + a(i,3) * v(3)
      enddo
      goto 99

!dir$ novector
 40   do i=1,nr
         w(i) = 
     +           a(i,1) * v(1)
     +         + a(i,2) * v(2)
     +         + a(i,3) * v(3)
     +         + a(i,4) * v(4)
      enddo
      goto 99

!dir$ novector
 50   do i=1,nr
         w(i) = 
     +           a(i,1) * v(1)
     +         + a(i,2) * v(2)
     +         + a(i,3) * v(3)
     +         + a(i,4) * v(4)
     +         + a(i,5) * v(5)
      enddo
      goto 99


!dir$ novector
 60   do i=1,nr
         w(i) = 
     +           a(i,1) * v(1)
     +         + a(i,2) * v(2)
     +         + a(i,3) * v(3)
     +         + a(i,4) * v(4)
     +         + a(i,5) * v(5)
     +         + a(i,6) * v(6)
      enddo
      goto 99

!dir$ novector
 70   do i=1,nr
         w(i) = 
     +           a(i,1) * v(1)
     +         + a(i,2) * v(2)
     +         + a(i,3) * v(3)
     +         + a(i,4) * v(4)
     +         + a(i,5) * v(5)
     +         + a(i,6) * v(6)
     +         + a(i,7) * v(7)
      enddo
      goto 99

!dir$ novector
 80   do i=1,nr
         w(i) = 
     +           a(i,1) * v(1)
     +         + a(i,2) * v(2)
     +         + a(i,3) * v(3)
     +         + a(i,4) * v(4)
     +         + a(i,5) * v(5)
     +         + a(i,6) * v(6)
     +         + a(i,7) * v(7)
     +         + a(i,8) * v(8)
      enddo
      irem=8
      goto 99

c --- finish up the remaining number of columns
 99   do j = irem+1,n,8
!dir$ novector
         do i = 1,nr
           w(i) =  w(i)
     +           + a(i,j  ) * v(j  )
     +           + a(i,j+1) * v(j+1)
     +           + a(i,j+2) * v(j+2)
     +           + a(i,j+3) * v(j+3)
     +           + a(i,j+4) * v(j+4)
     +           + a(i,j+5) * v(j+5)
     +           + a(i,j+6) * v(j+6)
     +           + a(i,j+7) * v(j+7)
         enddo
      enddo

c --- copy w array to output array v
!dir$ novector
      do k = 1,nr
         v(k) = w(k)
      enddo

      return
      end
