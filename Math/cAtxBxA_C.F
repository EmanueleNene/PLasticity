      subroutine cAtxBxA_C (a,b,c,na,nb,nc,n1,n2,work,con)
c --- primary function: do C = con*(At*B*A)
c
c  input arguments:
c     a        (dp,ar(na,*),in)    - matrix a
c     b        (dp,ar(nb,*),in)    - matrix b (must be square)
c     na       (int,sc,in)         - number of rows in matrix a
c     nb       (int,sc,in)         - number of rows in matrix b
c     nc       (int,sc,in)         - number of rows in matrix c
c     n1       (int,sc,in)         - number of rows in matrix c to fill
c     n2       (int,sc,in)         - number of columns in matrix c to fill
c     con      (dp,sc,in)          - multiplier on product added to sum
c
c  output arguments:
c     c        (dp,ar(nc,*),inout) - c = c + con*at*b*a
c     work     (dp,ar(n2,*),out)   - byproduct at*b
c
#include "impcom.inc"
c
      integer  na, nb, nc, n1, n2, i, j
      double precision   a(na,*),b(nb,*),c(nc,*), con, work(n2,*)

c --- second dimension of work must be the second dimension of matrix

c --- first do at * b, store result in work and store work as
c --- n2 by n1 matrix - that is, at by columns.

c --- if n2 > 64 eliminate all short loop directives
      if (n2.gt.64) goto 100

      goto (10,20,30,40,50,60,70,80,90), n1
        
  10  if (con .eq. 1.0d0) then
c --- multiply at * b - store in work by columns (n2 by n1)        
!dir$ ivdep,novector
            do j = 1,n2
               work(j,1) = 
     &                    + a(1,j)*b(1,1)
            enddo
      else
c --- multiply con*at * b - store in work by columns (n2 by n1)        
c --- con is a constant
!dir$ ivdep,novector
            do j = 1,n2
               work(j,1) = con*( 
     &                    + a(1,j)*b(1,1)
     &                                   )
            enddo
      endif

c --- multiply work * a store in c - lower triangle only
      do j = 1,n2
!dir$ ivdep,novector
         do i = j,n2
            c(i,j) =   c(i,j)
     &               + work(i,1)*a(1,j)
         enddo
      enddo
      goto 999 


 20   if (con .eq. 1.0d0) then
c --- multiply at * b - store in work by columns (n2 by n1)        
         do i = 1,2
!dir$ ivdep,novector
            do j = 1,n2
               work(j,i) = 
     &                    + a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
            enddo
         enddo
      else
c --- multiply con*at * b - store in work by columns (n2 by n1)        
c --- con is a constant
         do i = 1,2
!dir$ ivdep,novector
            do j = 1,n2
               work(j,i) = con*( 
     &                    + a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
     &                                   )
            enddo
         enddo
      endif

c --- multiply work * a store in c - lower triangle only
      do j = 1,n2
!dir$ ivdep,novector
         do i = j,n2
            c(i,j) = c(i,j)
     &               + work(i,1)*a(1,j)
     &               + work(i,2)*a(2,j)
         enddo
      enddo
      goto 999 


 30   if (con .eq. 1.0d0) then
c --- multiply at * b - store in work by columns (n2 by n1)        
         do i = 1,3
!dir$ ivdep,novector
            do j = 1,n2
               work(j,i) = 
     &                    + a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
     &                    + a(3,j)*b(3,i)
            enddo
         enddo
      else
c --- multiply con*at * b - store in work by columns (n2 by n1)        
c --- con is a constant
         do i = 1,3
!dir$ ivdep,novector
            do j = 1,n2
               work(j,i) = con*( 
     &                    + a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
     &                    + a(3,j)*b(3,i)
     &                                   )
            enddo
         enddo
      endif

c --- multiply work * a store in c - lower triangle only
      do j = 1,n2
!dir$ ivdep,novector
         do i = j,n2
            c(i,j) = c(i,j)
     &               + work(i,1)*a(1,j)
     &               + work(i,2)*a(2,j)
     &               + work(i,3)*a(3,j)
         enddo
      enddo
      goto 999 


 40   if (con .eq. 1.0d0) then
c --- multiply at * b - store in work by columns (n2 by n1)        
          do i = 1,4
!dir$ ivdep,novector
             do j = 1,n2
                work(j,i) = 
     &                     + a(1,j)*b(1,i)
     &                     + a(2,j)*b(2,i)
     &                     + a(3,j)*b(3,i)
     &                     + a(4,j)*b(4,i)
            enddo
         enddo
      else
c --- multiply con*at * b - store in work by columns (n2 by n1)        
c --- con is a constant
         do i = 1,4
!dir$ ivdep,novector
            do j = 1,n2
               work(j,i) = con*( 
     &                    + a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
     &                    + a(3,j)*b(3,i)
     &                    + a(4,j)*b(4,i)
     &                                   )
            enddo
         enddo
      endif

c --- multiply work * a store in c - lower triangle only
      do j = 1,n2
!dir$ ivdep,novector
         do i = j,n2
            c(i,j) = c(i,j)
     &               + work(i,1)*a(1,j)
     &               + work(i,2)*a(2,j)
     &               + work(i,3)*a(3,j)
     &               + work(i,4)*a(4,j)
         enddo
      enddo
      goto 999 


 50   if (con .eq. 1.0d0) then
c --- multiply at * b - store in work by columns (n2 by n1)        
          do i = 1,5
!dir$ ivdep,novector
             do j = 1,n2
                work(j,i) = 
     &                     + a(1,j)*b(1,i)
     &                     + a(2,j)*b(2,i)
     &                     + a(3,j)*b(3,i)
     &                     + a(4,j)*b(4,i)
     &                     + a(5,j)*b(5,i)
             enddo
          enddo
      else
c --- multiply con*at * b - store in work by columns (n2 by n1)        
c --- con is a constant
         do i = 1,5
            do j = 1,n2
               work(j,i) = con*( 
     &                    + a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
     &                    + a(3,j)*b(3,i)
     &                    + a(4,j)*b(4,i)
     &                    + a(5,j)*b(5,i)
     &                                   )
            enddo
         enddo
      endif

c --- multiply work * a store in c - lower triangle only
      do j = 1,n2
!dir$ ivdep,novector
         do i = j,n2
            c(i,j) =   c(i,j)
     &               + work(i,1)*a(1,j)
     &               + work(i,2)*a(2,j)
     &               + work(i,3)*a(3,j)
     &               + work(i,4)*a(4,j)
     &               + work(i,5)*a(5,j)
         enddo
      enddo
      goto 999 


 60   if (con .eq. 1.0d0) then
c --- multiply at * b - store in work by columns (n2 by n1)        
         do i = 1,6
!dir$ ivdep,novector
            do j = 1,n2
               work(j,i) = 
     &                    + a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
     &                    + a(3,j)*b(3,i)
     &                    + a(4,j)*b(4,i)
     &                    + a(5,j)*b(5,i)
     &                    + a(6,j)*b(6,i)
            enddo
         enddo
      else
c --- multiply con*at * b - store in work by columns (n2 by n1)        
c --- con is a constant
         do i = 1,6
!dir$ ivdep,novector
            do j = 1,n2
               work(j,i) = con*( 
     &                    + a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
     &                    + a(3,j)*b(3,i)
     &                    + a(4,j)*b(4,i)
     &                    + a(5,j)*b(5,i)
     &                    + a(6,j)*b(6,i)
     &                                   )
            enddo
         enddo
      endif

c --- multiply work * a store in c - lower triangle only
      do j = 1,n2
!dir$ ivdep,novector
         do i = j,n2
            c(i,j) =   c(i,j)
     &               + work(i,1)*a(1,j)
     &               + work(i,2)*a(2,j)
     &               + work(i,3)*a(3,j)
     &               + work(i,4)*a(4,j)
     &               + work(i,5)*a(5,j)
     &               + work(i,6)*a(6,j)
         enddo
      enddo
      goto 999 


 70   if (con .eq. 1.0d0) then
c --- multiply at * b - store in work by columns (n2 by n1)        
         do i = 1,7
!dir$ ivdep,novector
            do j = 1,n2
               work(j,i) = 
     &                    + a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
     &                    + a(3,j)*b(3,i)
     &                    + a(4,j)*b(4,i)
     &                    + a(5,j)*b(5,i)
     &                    + a(6,j)*b(6,i)
     &                    + a(7,j)*b(7,i)
            enddo
         enddo
      else
c --- multiply con*at * b - store in work by columns (n2 by n1)        
c --- con is a constant
         do i = 1,7
!dir$ ivdep,novector
            do j = 1,n2
               work(j,i) = con*( 
     &                    + a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
     &                    + a(3,j)*b(3,i)
     &                    + a(4,j)*b(4,i)
     &                    + a(5,j)*b(5,i)
     &                    + a(6,j)*b(6,i)
     &                    + a(7,j)*b(7,i)
     &                                   )
            enddo
         enddo
      endif

c --- multiply work * a store in c - lower triangle only
      do j = 1,n2
!dir$ ivdep,novector
         do i = j,n2
            c(i,j) = c(i,j)
     &               + work(i,1)*a(1,j)
     &               + work(i,2)*a(2,j)
     &               + work(i,3)*a(3,j)
     &               + work(i,4)*a(4,j)
     &               + work(i,5)*a(5,j)
     &               + work(i,6)*a(6,j)
     &               + work(i,7)*a(7,j)
         enddo
      enddo
      goto 999 


 80   if (con .eq. 1.0d0) then
c --- multiply at * b - store in work by columns (n2 by n1)        
         do i = 1,8
!dir$ ivdep,novector
            do j = 1,n2
               work(j,i) = 
     &                    + a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
     &                    + a(3,j)*b(3,i)
     &                    + a(4,j)*b(4,i)
     &                    + a(5,j)*b(5,i)
     &                    + a(6,j)*b(6,i)
     &                    + a(7,j)*b(7,i)
     &                    + a(8,j)*b(8,i)
            enddo
         enddo
      else
c --- multiply con*at * b - store in work by columns (n2 by n1)        
c --- con is a constant
         do i = 1,8
!dir$ ivdep,novector
            do j = 1,n2
               work(j,i) = con*( 
     &                    + a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
     &                    + a(3,j)*b(3,i)
     &                    + a(4,j)*b(4,i)
     &                    + a(5,j)*b(5,i)
     &                    + a(6,j)*b(6,i)
     &                    + a(7,j)*b(7,i)
     &                    + a(8,j)*b(8,i)
     &                                   )
            enddo
         enddo
      endif

c --- multiply work * a store in c - lower triangle only
      do j = 1,n2
!dir$ ivdep,novector
         do i = j,n2
            c(i,j) =   c(i,j)
     &               + work(i,1)*a(1,j)
     &               + work(i,2)*a(2,j)
     &               + work(i,3)*a(3,j)
     &               + work(i,4)*a(4,j)
     &               + work(i,5)*a(5,j)
     &               + work(i,6)*a(6,j)
     &               + work(i,7)*a(7,j)
     &               + work(i,8)*a(8,j)
         enddo
      enddo
      goto 999 


 90   if (con.eq.1) then
c --- multiply at * b - store in work by columns (n2 by n1)        
         do i = 1,9
!dir$ ivdep,novector
            do j = 1,n2
               work(j,i) = 
     &                    + a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
     &                    + a(3,j)*b(3,i)
     &                    + a(4,j)*b(4,i)
     &                    + a(5,j)*b(5,i)
     &                    + a(6,j)*b(6,i)
     &                    + a(7,j)*b(7,i)
     &                    + a(8,j)*b(8,i)
     &                    + a(9,j)*b(9,i)
            enddo
         enddo
      else
c --- multiply con*at * b - store in work by columns (n2 by n1)        
c --- con is a constant
         do i = 1,9
!dir$ ivdep,novector
            do j = 1,n2
               work(j,i) = con*( 
     &                    + a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
     &                    + a(3,j)*b(3,i)
     &                    + a(4,j)*b(4,i)
     &                    + a(5,j)*b(5,i)
     &                    + a(6,j)*b(6,i)
     &                    + a(7,j)*b(7,i)
     &                    + a(8,j)*b(8,i)
     &                    + a(9,j)*b(9,i)
     &                                   )
            enddo
         enddo
      endif

c --- multiply work * a store in c - lower triangle only
      do j = 1,n2
!dir$ ivdep,novector
         do i = j,n2
            c(i,j) = c(i,j)
     &               + work(i,1)*a(1,j)
     &               + work(i,2)*a(2,j)
     &               + work(i,3)*a(3,j)
     &               + work(i,4)*a(4,j)
     &               + work(i,5)*a(5,j)
     &               + work(i,6)*a(6,j)
     &               + work(i,7)*a(7,j)
     &               + work(i,8)*a(8,j)
     &               + work(i,9)*a(9,j)
         enddo
      enddo
      goto 999 


 100  goto (110,120,130,140,150,160,170,180,190) n1
        
 110  if (con .eq. 1.0d0) then
c --- multiply at * b - store in work by columns (n2 by n1)        
            do j = 1,n2
               work(j,1) = a(1,j)*b(1,1)
            enddo
      else
c --- multiply con*at * b - store in work by columns (n2 by n1)        
c --- con is a constant
            do j = 1,n2
               work(j,1) = con*( 
     &                     + a(1,j)*b(1,1))
            enddo
      endif

c --- multiply work * a store in c - lower triangle only
      do j = 1,n2
         do i = j,n2
            c(i,j) = c(i,j)
     &               + work(i,1)*a(1,j)
         enddo
      enddo
      goto 999 


 120  if (con .eq. 1.0d0) then
c --- multiply at * b - store in work by columns (n2 by n1)        
         do i = 1,2
            do j = 1,n2
               work(j,i) =  a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
            enddo
         enddo
      else
c --- multiply con*at * b - store in work by columns (n2 by n1)        
c --- con is a constant
         do i = 1,2
            do j = 1,n2
               work(j,i) = con*( 
     &                    + a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
     &                                   )
            enddo
         enddo
      endif

c --- multiply work * a store in c - lower triangle only
      do j = 1,n2
         do i = j,n2
            c(i,j) = c(i,j)
     &               + work(i,1)*a(1,j)
     &               + work(i,2)*a(2,j)
         enddo
      enddo
      goto 999 


 130  if (con .eq. 1.0d0) then
c --- multiply at * b - store in work by columns (n2 by n1)        
         do i = 1,3
            do j = 1,n2
               work(j,i) =  a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
     &                    + a(3,j)*b(3,i)
            enddo
         enddo
      else
c --- multiply con*at * b - store in work by columns (n2 by n1)        
c --- con is a constant
         do i = 1,3
            do j = 1,n2
               work(j,i) = con*( 
     &                    + a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
     &                    + a(3,j)*b(3,i)
     &                                   )
            enddo
         enddo
      endif

c --- multiply work * a store in c - lower triangle only
      do j = 1,n2
         do i = j,n2
            c(i,j) = c(i,j)
     &               + work(i,1)*a(1,j)
     &               + work(i,2)*a(2,j)
     &               + work(i,3)*a(3,j)
         enddo
      enddo
      goto 999 


 140  if (con .eq. 1.0d0) then
c --- multiply at * b - store in work by columns (n2 by n1)        
         do i = 1,4
            do j = 1,n2
               work(j,i) =  a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
     &                    + a(3,j)*b(3,i)
     &                    + a(4,j)*b(4,i)
            enddo
         enddo
      else
c --- multiply con*at * b - store in work by columns (n2 by n1)        
c --- con is a constant
         do i = 1,4
            do j = 1,n2
               work(j,i) = con*( 
     &                    + a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
     &                    + a(3,j)*b(3,i)
     &                    + a(4,j)*b(4,i)
     &                                   )
            enddo
         enddo
      endif

c --- multiply work * a store in c - lower triangle only
      do j = 1,n2
         do i = j,n2
            c(i,j) = c(i,j)
     &               + work(i,1)*a(1,j)
     &               + work(i,2)*a(2,j)
     &               + work(i,3)*a(3,j)
     &               + work(i,4)*a(4,j)
         enddo
      enddo
      goto 999 


 150  if (con .eq. 1.0d0) then
c --- multiply at * b - store in work by columns (n2 by n1)        
         do i = 1,5
            do j = 1,n2
               work(j,i) =  a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
     &                    + a(3,j)*b(3,i)
     &                    + a(4,j)*b(4,i)
     &                    + a(5,j)*b(5,i)
            enddo
         enddo
      else
c --- multiply con*at * b - store in work by columns (n2 by n1)        
c --- con is a constant
         do i = 1,5
            do j = 1,n2
               work(j,i) = con*( 
     &                    + a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
     &                    + a(3,j)*b(3,i)
     &                    + a(4,j)*b(4,i)
     &                    + a(5,j)*b(5,i)
     &                                   )
            enddo
         enddo
      endif

c --- multiply work * a store in c - lower triangle only
      do j = 1,n2
         do i = j,n2
            c(i,j) = c(i,j)
     &               + work(i,1)*a(1,j)
     &               + work(i,2)*a(2,j)
     &               + work(i,3)*a(3,j)
     &               + work(i,4)*a(4,j)
     &               + work(i,5)*a(5,j)
         enddo
      enddo
      goto 999 


 160  if (con .eq. 1.0d0) then
c --- multiply at * b - store in work by columns (n2 by n1)        
         do i = 1,6
            do j = 1,n2
               work(j,i) =  a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
     &                    + a(3,j)*b(3,i)
     &                    + a(4,j)*b(4,i)
     &                    + a(5,j)*b(5,i)
     &                    + a(6,j)*b(6,i)
            enddo
         enddo
      else
c --- multiply con*at * b - store in work by columns (n2 by n1)        
c --- con is a constant
         do i = 1,6
            do j = 1,n2
               work(j,i) = con*( 
     &                    + a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
     &                    + a(3,j)*b(3,i)
     &                    + a(4,j)*b(4,i)
     &                    + a(5,j)*b(5,i)
     &                    + a(6,j)*b(6,i)
     &                                   )
            enddo
         enddo
      endif

c --- multiply work * a store in c - lower triangle only
      do j = 1,n2
         do i = j,n2
            c(i,j) = c(i,j)
     &               + work(i,1)*a(1,j)
     &               + work(i,2)*a(2,j)
     &               + work(i,3)*a(3,j)
     &               + work(i,4)*a(4,j)
     &               + work(i,5)*a(5,j)
     &               + work(i,6)*a(6,j)
         enddo
      enddo
      goto 999 


 170  if (con .eq. 1.0d0) then
c --- multiply at * b - store in work by columns (n2 by n1)        
         do i = 1,7
            do j = 1,n2
               work(j,i) =  a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
     &                    + a(3,j)*b(3,i)
     &                    + a(4,j)*b(4,i)
     &                    + a(5,j)*b(5,i)
     &                    + a(6,j)*b(6,i)
     &                    + a(7,j)*b(7,i)
            enddo
         enddo
      else
c --- multiply con*at * b - store in work by columns (n2 by n1)        
c --- con is a constant
         do i = 1,7
            do j = 1,n2
               work(j,i) = con*( 
     &                    + a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
     &                    + a(3,j)*b(3,i)
     &                    + a(4,j)*b(4,i)
     &                    + a(5,j)*b(5,i)
     &                    + a(6,j)*b(6,i)
     &                    + a(7,j)*b(7,i)
     &                                   )
            enddo
         enddo
      endif

c --- multiply work * a store in c - lower triangle only
      do j = 1,n2
         do i = j,n2
            c(i,j) = c(i,j)
     &               + work(i,1)*a(1,j)
     &               + work(i,2)*a(2,j)
     &               + work(i,3)*a(3,j)
     &               + work(i,4)*a(4,j)
     &               + work(i,5)*a(5,j)
     &               + work(i,6)*a(6,j)
     &               + work(i,7)*a(7,j)
         enddo
      enddo
      goto 999 


 180  if (con .eq. 1.0d0) then
c --- multiply at * b - store in work by columns (n2 by n1)        
         do i = 1,8
            do j = 1,n2
               work(j,i) =  a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
     &                    + a(3,j)*b(3,i)
     &                    + a(4,j)*b(4,i)
     &                    + a(5,j)*b(5,i)
     &                    + a(6,j)*b(6,i)
     &                    + a(7,j)*b(7,i)
     &                    + a(8,j)*b(8,i)
            enddo
         enddo
      else
c --- multiply con*at * b - store in work by columns (n2 by n1)        
c --- con is a constant
         do i = 1,8
            do j = 1,n2
               work(j,i) = con*( 
     &                    + a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
     &                    + a(3,j)*b(3,i)
     &                    + a(4,j)*b(4,i)
     &                    + a(5,j)*b(5,i)
     &                    + a(6,j)*b(6,i)
     &                    + a(7,j)*b(7,i)
     &                    + a(8,j)*b(8,i)
     &                                   )
            enddo
         enddo
      endif

c --- multiply work * a store in c - lower triangle only
      do j = 1,n2
         do i = j,n2
            c(i,j) = c(i,j)
     &               + work(i,1)*a(1,j)
     &               + work(i,2)*a(2,j)
     &               + work(i,3)*a(3,j)
     &               + work(i,4)*a(4,j)
     &               + work(i,5)*a(5,j)
     &               + work(i,6)*a(6,j)
     &               + work(i,7)*a(7,j)
     &               + work(i,8)*a(8,j)
         enddo
      enddo
      goto 999 


 190  if (con .eq. 1.0d0) then
c --- multiply at * b - store in work by columns (n2 by n1)        
         do i = 1,9
            do j = 1,n2
               work(j,i) =  a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
     &                    + a(3,j)*b(3,i)
     &                    + a(4,j)*b(4,i)
     &                    + a(5,j)*b(5,i)
     &                    + a(6,j)*b(6,i)
     &                    + a(7,j)*b(7,i)
     &                    + a(8,j)*b(8,i)
     &                    + a(9,j)*b(9,i)
            enddo
         enddo
      else
c --- multiply con*at * b - store in work by columns (n2 by n1)        
c --- con is a constant
         do i = 1,9
            do j = 1,n2
               work(j,i) = con*( 
     &                    + a(1,j)*b(1,i)
     &                    + a(2,j)*b(2,i)
     &                    + a(3,j)*b(3,i)
     &                    + a(4,j)*b(4,i)
     &                    + a(5,j)*b(5,i)
     &                    + a(6,j)*b(6,i)
     &                    + a(7,j)*b(7,i)
     &                    + a(8,j)*b(8,i)
     &                    + a(9,j)*b(9,i)
     &                                   )
            enddo
         enddo
      endif

c --- multiply work * a store in c - lower triangle only
      do j = 1,n2
         do i = j,n2
            c(i,j) = c(i,j)
     &               + work(i,1)*a(1,j)
     &               + work(i,2)*a(2,j)
     &               + work(i,3)*a(3,j)
     &               + work(i,4)*a(4,j)
     &               + work(i,5)*a(5,j)
     &               + work(i,6)*a(6,j)
     &               + work(i,7)*a(7,j)
     &               + work(i,8)*a(8,j)
     &               + work(i,9)*a(9,j)
         enddo
      enddo

 999  return
      end
