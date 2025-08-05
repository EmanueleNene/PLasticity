      subroutine BtxAxB_A (a,na,b,nb,nold,nnew)
c --- primary function: do  A = Bt * A * B ,  where A is symmetric
c
c  input arguments:
c     na       (int,sc,in)            - first dimension of the a matrix
c     b        (dp,ar(nb,nnew),in)    - matrix to pre and post multiply a by
c                                       (part used may be rectangular(nold x nnew))
c     nb       (int,sc,in)            - first dimension of the c matrix
c     nold     (int,sc,in)            - size of part of 'A' matrix that is to be processed(input size)
c     nnew     (int,sc,in)            - size of part of 'A' matrix that results from this operation(output size)
c
c  output arguments:
c     a        (dp,ar(na,na),inout)   - matrix to be pre and post multiplied
c                                       part operated on must be square(nold x nold) and symmetric
c                                       resulting matrix (still square(nnew x nnew) and symmetric)
c
#include "impcom.inc"
c
      integer  na, nb, nold, nnew,  i, j
      double precision  a(na,na), b(nb,nnew)
      integer  irem, ist, k
      double precision  w(64,64)

c --- first do ct x a - store in array w by columns
      irem = mod(nold,8)

      goto (80,10,20,30,40,50,60,70) irem+1


 10   do i=1,nold
!dir$ ivdep,novector
         do j=1,nnew
            w(j,i) =
     +                    b(1,j) * a(1,i) 
         enddo
      enddo
      ist=irem+1
      goto 90

 20   do i=1,nold
!dir$ ivdep,novector
         do j=1,nnew
            w(j,i) =
     +                    b(1,j) * a(1,i) 
     +                  + b(2,j) * a(2,i) 
         enddo
      enddo
      ist=irem+1
      goto 90

 30   do i=1,nold
!dir$ ivdep,novector
         do j=1,nnew
            w(j,i) =
     +                    b(1,j) * a(1,i) 
     +                  + b(2,j) * a(2,i) 
     +                  + b(3,j) * a(3,i) 
         enddo
      enddo
      ist=irem+1
      goto 90

 40   do i=1,nold
!dir$ ivdep,novector
         do j=1,nnew
            w(j,i) =
     +                    b(1,j) * a(1,i) 
     +                  + b(2,j) * a(2,i) 
     +                  + b(3,j) * a(3,i) 
     +                  + b(4,j) * a(4,i) 
         enddo
      enddo
      ist=irem+1
      goto 90

 50   do i=1,nold
!dir$ ivdep,novector
         do j=1,nnew
            w(j,i) =
     +                    b(1,j) * a(1,i) 
     +                  + b(2,j) * a(2,i) 
     +                  + b(3,j) * a(3,i) 
     +                  + b(4,j) * a(4,i) 
     +                  + b(5,j) * a(5,i) 
         enddo
      enddo
      ist=irem+1
      goto 90

 60   do i=1,nold
!dir$ ivdep,novector
         do j=1,nnew
            w(j,i) =
     +                    b(1,j) * a(1,i) 
     +                  + b(2,j) * a(2,i) 
     +                  + b(3,j) * a(3,i) 
     +                  + b(4,j) * a(4,i) 
     +                  + b(5,j) * a(5,i) 
     +                  + b(6,j) * a(6,i) 
         enddo
      enddo
      ist=irem+1
      goto 90

 70   do i=1,nold
!dir$ ivdep,novector
         do j=1,nnew
            w(j,i) =
     +                    b(1,j) * a(1,i) 
     +                  + b(2,j) * a(2,i) 
     +                  + b(3,j) * a(3,i) 
     +                  + b(4,j) * a(4,i) 
     +                  + b(5,j) * a(5,i) 
     +                  + b(6,j) * a(6,i) 
     +                  + b(7,j) * a(7,i) 
         enddo
      enddo
      ist=irem+1
      goto 90

 80   do i=1,nold
!dir$ ivdep,novector
         do j=1,nnew
            w(j,i) =
     +                    b(1,j) * a(1,i) 
     +                  + b(2,j) * a(2,i) 
     +                  + b(3,j) * a(3,i) 
     +                  + b(4,j) * a(4,i) 
     +                  + b(5,j) * a(5,i) 
     +                  + b(6,j) * a(6,i) 
     +                  + b(7,j) * a(7,i) 
     +                  + b(8,j) * a(8,i) 
         enddo
      enddo
      ist=9
      
 90   do k=ist,nold,8
         do i=1,nold
!dir$ ivdep,novector
            do j=1,nnew
                 w(j,i) = w(j,i)
     +                      + b(k  ,j) * a(k  ,i)
     +                      + b(k+1,j) * a(k+1,i)
     +                      + b(k+2,j) * a(k+2,i)
     +                      + b(k+3,j) * a(k+3,i)
     +                      + b(k+4,j) * a(k+4,i)
     +                      + b(k+5,j) * a(k+5,i)
     +                      + b(k+6,j) * a(k+6,i)
     +                      + b(k+7,j) * a(k+7,i)
            enddo
         enddo
      enddo

c --- now finish with w x c
c --- for larger matrices compute only half of the result and
c --- then copy the symmetric half - for smaller matrices it is
c --- faster to compute in one piece

      if ((nnew.gt.22).and.(nold.gt.8)) goto 200
      goto (180,110,120,130,140,150,160,170) irem+1

 110  do j=1,nnew
!dir$ ivdep,novector
         do i=1,nnew
            a(i,j) =
     +                  w(i,1) * b(1,j)
         enddo
      enddo
      goto 190
 
 120  do j=1,nnew
!dir$ ivdep,novector
         do i=1,nnew
            a(i,j) =
     +                  w(i,1) * b(1,j)
     +                + w(i,2) * b(2,j)
         enddo
      enddo
      goto 190

 130  do j=1,nnew
!dir$ ivdep,novector
         do i=1,nnew
            a(i,j) =
     +                  w(i,1) * b(1,j)
     +                + w(i,2) * b(2,j)
     +                + w(i,3) * b(3,j)
         enddo
      enddo
      goto 190


 140  do j=1,nnew
!dir$ ivdep,novector
         do i=1,nnew
            a(i,j) =
     +                  w(i,1) * b(1,j)
     +                + w(i,2) * b(2,j)
     +                + w(i,3) * b(3,j)
     +                + w(i,4) * b(4,j)
         enddo
      enddo
      goto 190


 150  do j=1,nnew
!dir$ ivdep,novector
         do i=1,nnew
            a(i,j) =
     +                  w(i,1) * b(1,j)
     +                + w(i,2) * b(2,j)
     +                + w(i,3) * b(3,j)
     +                + w(i,4) * b(4,j)
     +                + w(i,5) * b(5,j)
         enddo
      enddo
      goto 190

 160  do j=1,nnew
!dir$ ivdep,novector
         do i=1,nnew
            a(i,j) =
     +                  w(i,1) * b(1,j)
     +                + w(i,2) * b(2,j)
     +                + w(i,3) * b(3,j)
     +                + w(i,4) * b(4,j)
     +                + w(i,5) * b(5,j)
     +                + w(i,6) * b(6,j)
         enddo
      enddo
      goto 190

 170  do j=1,nnew
!dir$ ivdep,novector
         do i=1,nnew
            a(i,j) =
     +                  w(i,1) * b(1,j)
     +                + w(i,2) * b(2,j)
     +                + w(i,3) * b(3,j)
     +                + w(i,4) * b(4,j)
     +                + w(i,5) * b(5,j)
     +                + w(i,6) * b(6,j)
     +                + w(i,7) * b(7,j)
         enddo
      enddo
      goto 190

 180  do j=1,nnew
!dir$ ivdep,novector
         do i=1,nnew
            a(i,j) = 
     +                  w(i,1) * b(1,j)
     +                + w(i,2) * b(2,j)
     +                + w(i,3) * b(3,j)
     +                + w(i,4) * b(4,j)
     +                + w(i,5) * b(5,j)
     +                + w(i,6) * b(6,j)
     +                + w(i,7) * b(7,j)
     +                + w(i,8) * b(8,j)
         enddo
      enddo

 190  do k=ist,nold,8
         do j=1,nnew
!dir$ ivdep,novector
            do i=1,nnew
               a(i,j) = a(i,j)
     +                   + w(i,k  ) * b(k  ,j)
     +                   + w(i,k+1) * b(k+1,j)
     +                   + w(i,k+2) * b(k+2,j)
     +                   + w(i,k+3) * b(k+3,j)
     +                   + w(i,k+4) * b(k+4,j)
     +                   + w(i,k+5) * b(k+5,j)
     +                   + w(i,k+6) * b(k+6,j)
     +                   + w(i,k+7) * b(k+7,j)
            enddo
         enddo
      enddo
      goto 999 


c --- this section computes only half of the symmetric output matrix
c --- and then fills in the rest at the end - used only for larger
c --- matrices

 200  goto (280,210,220,230,240,250,260,270) irem+1

 210  do j=1,nnew
!dir$ ivdep,novector
         do i=j,nnew
            a(i,j) =
     +                  w(i,1) * b(1,j)
         enddo
      enddo
      goto 290
 
 220  do j=1,nnew
!dir$ ivdep,novector
         do i=j,nnew
            a(i,j) =
     +                  w(i,1) * b(1,j)
     +                + w(i,2) * b(2,j)
         enddo
      enddo
      goto 290

 230  do j=1,nnew
!dir$ ivdep,novector
         do i=j,nnew
            a(i,j) =
     +                  w(i,1) * b(1,j)
     +                + w(i,2) * b(2,j)
     +                + w(i,3) * b(3,j)
         enddo
      enddo
      goto 290

 240  do j=1,nnew
!dir$ ivdep,novector
         do i=j,nnew
            a(i,j) =
     +                  w(i,1) * b(1,j)
     +                + w(i,2) * b(2,j)
     +                + w(i,3) * b(3,j)
     +                + w(i,4) * b(4,j)
         enddo
      enddo
      goto 290


 250  do j=1,nnew
!dir$ ivdep,novector
         do i=j,nnew
            a(i,j) =
     +                  w(i,1) * b(1,j)
     +                + w(i,2) * b(2,j)
     +                + w(i,3) * b(3,j)
     +                + w(i,4) * b(4,j)
     +                + w(i,5) * b(5,j)
         enddo
      enddo
      goto 290

 260  do j=1,nnew
!dir$ ivdep,novector
         do i=j,nnew
            a(i,j) =
     +                  w(i,1) * b(1,j)
     +                + w(i,2) * b(2,j)
     +                + w(i,3) * b(3,j)
     +                + w(i,4) * b(4,j)
     +                + w(i,5) * b(5,j)
     +                + w(i,6) * b(6,j)
         enddo
      enddo
      goto 290

 270  do j=1,nnew
!dir$ ivdep,novector
         do i=j,nnew
            a(i,j) =
     +                  w(i,1) * b(1,j)
     +                + w(i,2) * b(2,j)
     +                + w(i,3) * b(3,j)
     +                + w(i,4) * b(4,j)
     +                + w(i,5) * b(5,j)
     +                + w(i,6) * b(6,j)
     +                + w(i,7) * b(7,j)
         enddo
      enddo
      goto 290

 280  do j=1,nnew
!dir$ ivdep,novector
         do i=j,nnew
            a(i,j) = 
     +                  w(i,1) * b(1,j)
     +                + w(i,2) * b(2,j)
     +                + w(i,3) * b(3,j)
     +                + w(i,4) * b(4,j)
     +                + w(i,5) * b(5,j)
     +                + w(i,6) * b(6,j)
     +                + w(i,7) * b(7,j)
     +                + w(i,8) * b(8,j)
         enddo
      enddo

 290  do k=ist,nold,8
         do j=1,nnew
!dir$ ivdep,novector
            do i=j,nnew
               a(i,j) = a(i,j)
     +                   + w(i,k  ) * b(k  ,j)
     +                   + w(i,k+1) * b(k+1,j)
     +                   + w(i,k+2) * b(k+2,j)
     +                   + w(i,k+3) * b(k+3,j)
     +                   + w(i,k+4) * b(k+4,j)
     +                   + w(i,k+5) * b(k+5,j)
     +                   + w(i,k+6) * b(k+6,j)
     +                   + w(i,k+7) * b(k+7,j)
            enddo
         enddo
      enddo

      do j=2,nnew
!dir$ ivdep,novector
         do i=1,j-1
            a(i,j)=a(j,i)
         enddo
      enddo

 999  return
      end
