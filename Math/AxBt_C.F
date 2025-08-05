      subroutine AxBt_C (a,b,c, na,nb,nc, n1,n2,n3)
c --- primary function:  (a) * (b)t = (c)     t means transpose
c                        summs in c if n3.gt.0
c
c  input arguments:
c     a        (dp,ar(na,*),in)   - matrix a
c     b        (dp,ar(nb,*),in)   - matrix b
c     na       (int,sc,in)        - number of rows in matrix a
c     nb       (int,sc,in)        - number of rows in matrix b
c     nc       (int,sc,in)        - number of rows in matrix c
c     n1       (int,sc,in)        - number of rows in matrix c to fill
c     n2       (int,sc,in)        - number of columns in matrix c to fill
c     n3       (int,sc,in)        - number of columns in matrix a and
c                                   number of rows in matrix b
c                                   to work with (the two need
c                                   to be the same for the inner product)
c
c  output arguments:
c     c        (dp,ar(nc,*),out)  - product matrix c
c
#include "impcom.inc"
c
      integer  na,nb,nc, n1,n2,n3, n, i,j,k, irem
      double precision   a(na,*),b(nb,*),c(nc,*)
      integer  j1, k1, k2, k3, k4, k5, k6, k7

      n = - n3
      if (n .le. 0) then
c --- new target matrix

         n = n3
         irem = mod(n,8)
         goto (80,10,20,30,40,50,60,70) irem + 1
        
 10      do j = 1,n2
            do i = 1,n1
               c(i,j) =   a(i,1)*b(j,1)
            enddo
         enddo
         goto 99

 20      do j = 1,n2
            do i = 1,n1
               c(i,j) =   a(i,1)*b(j,1)
     &                  + a(i,2)*b(j,2)
            enddo
         enddo
         goto 99

 30      do j = 1,n2
            do i = 1,n1
               c(i,j) =   a(i,1)*b(j,1)
     &                  + a(i,2)*b(j,2)
     &                  + a(i,3)*b(j,3)
            enddo
         enddo
         goto 99

 40      do j = 1,n2
            do i = 1,n1
               c(i,j) =   a(i,1)*b(j,1)
     &                  + a(i,2)*b(j,2)
     &                  + a(i,3)*b(j,3)
     &                  + a(i,4)*b(j,4)
            enddo
         enddo
         goto 99

 50      do j = 1,n2
            do i = 1,n1
               c(i,j) =   a(i,1)*b(j,1)
     &                  + a(i,2)*b(j,2)
     &                  + a(i,3)*b(j,3)
     &                  + a(i,4)*b(j,4)
     &                  + a(i,5)*b(j,5)
            enddo
         enddo
         goto 99

 60      do j = 1,n2
            do i = 1,n1
               c(i,j) =   a(i,1)*b(j,1)
     &                  + a(i,2)*b(j,2)
     &                  + a(i,3)*b(j,3)
     &                  + a(i,4)*b(j,4)
     &                  + a(i,5)*b(j,5)
     &                  + a(i,6)*b(j,6)
            enddo
         enddo
         goto 99

 70      do j = 1,n2
            do i = 1,n1
               c(i,j) =   a(i,1)*b(j,1)
     &                  + a(i,2)*b(j,2)
     &                  + a(i,3)*b(j,3)
     &                  + a(i,4)*b(j,4)
     &                  + a(i,5)*b(j,5)
     &                  + a(i,6)*b(j,6)
     &                  + a(i,7)*b(j,7)
            enddo
         enddo
         goto 99

 80      irem = 8
         do j = 1,n2
            do i = 1,n1
               c(i,j) =   a(i,1)*b(j,1)
     &                  + a(i,2)*b(j,2)
     &                  + a(i,3)*b(j,3)
     &                  + a(i,4)*b(j,4)
     &                  + a(i,5)*b(j,5)
     &                  + a(i,6)*b(j,6)
     &                  + a(i,7)*b(j,7)
     &                  + a(i,8)*b(j,8)
            enddo
         enddo
         goto 99
         
      else


         irem = mod(n,8)
         goto (99,110,120,130,140,150,160,170) irem + 1

 110     do j = 1,n2
            do i = 1,n1
               c(i,j) = c(i,j) + a(i,1)*b(j,1)
            enddo
         enddo
         goto 99

 120     do j = 1,n2
            do i = 1,n1
               c(i,j) = c(i,j) + a(i,1)*b(j,1)
     &                         + a(i,2)*b(j,2)
            enddo
         enddo
         goto 99

 130     do j = 1,n2
            do i = 1,n1
               c(i,j) = c(i,j) + a(i,1)*b(j,1)
     &                         + a(i,2)*b(j,2)
     &                         + a(i,3)*b(j,3)
            enddo
         enddo
         goto 99

 140     do j = 1,n2
            do i = 1,n1
               c(i,j) = c(i,j) + a(i,1)*b(j,1)
     &                         + a(i,2)*b(j,2)
     &                         + a(i,3)*b(j,3)
     &                         + a(i,4)*b(j,4)
            enddo
         enddo
         goto 99


 150     do j = 1,n2
            do i = 1,n1
               c(i,j) = c(i,j) + a(i,1)*b(j,1)
     &                         + a(i,2)*b(j,2)
     &                         + a(i,3)*b(j,3)
     &                         + a(i,4)*b(j,4)
     &                         + a(i,5)*b(j,5)
            enddo
         enddo
         goto 99

 160     do j = 1,n2
            do i = 1,n1
               c(i,j) = c(i,j) + a(i,1)*b(j,1)
     &                         + a(i,2)*b(j,2)
     &                         + a(i,3)*b(j,3)
     &                         + a(i,4)*b(j,4)
     &                         + a(i,5)*b(j,5)
     &                         + a(i,6)*b(j,6)
            enddo
         enddo
         goto 99

 170     do j = 1,n2
            do i = 1,n1
               c(i,j) = c(i,j) + a(i,1)*b(j,1)
     &                         + a(i,2)*b(j,2)
     &                         + a(i,3)*b(j,3)
     &                         + a(i,4)*b(j,4)
     &                         + a(i,5)*b(j,5)
     &                         + a(i,6)*b(j,6)
     &                         + a(i,7)*b(j,7)
            enddo
         enddo

      endif

 99   do k = irem+1,n,8
         k1 = k + 1
         k2 = k + 2
         k3 = k + 3
         k4 = k + 4
         k5 = k + 5
         k6 = k + 6
         k7 = k + 7
         do j = 1,n2-1,2
            j1 = j + 1
!dir$ ivdep
            do i = 1,n1
               c(i,j) = c(i,j) + a(i,k )*b(j,k )
     &                         + a(i,k1)*b(j,k1)
     &                         + a(i,k2)*b(j,k2)
     &                         + a(i,k3)*b(j,k3)
     &                         + a(i,k4)*b(j,k4)
     &                         + a(i,k5)*b(j,k5)
     &                         + a(i,k6)*b(j,k6)
     &                         + a(i,k7)*b(j,k7)
               c(i,j1) = c(i,j1) + a(i,k )*b(j1,k )
     &                           + a(i,k1)*b(j1,k1)
     &                           + a(i,k2)*b(j1,k2)
     &                           + a(i,k3)*b(j1,k3)
     &                           + a(i,k4)*b(j1,k4)
     &                           + a(i,k5)*b(j1,k5)
     &                           + a(i,k6)*b(j1,k6)
     &                           + a(i,k7)*b(j1,k7)
            enddo
         enddo
c
         do j = j,n2
            do i = 1,n1
               c(i,j) = c(i,j) + a(i,k )*b(j,k )
     &                         + a(i,k1)*b(j,k1)
     &                         + a(i,k2)*b(j,k2)
     &                         + a(i,k3)*b(j,k3)
     &                         + a(i,k4)*b(j,k4)
     &                         + a(i,k5)*b(j,k5)
     &                         + a(i,k6)*b(j,k6)
     &                         + a(i,k7)*b(j,k7)
            enddo
         enddo
      enddo

      return
      end
