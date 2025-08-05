      subroutine CalCmt(cmtr,feba,deba,nnod,nr)
c --- primary function: compute c matrix for stress stiffnessing
c
c     ref: nour-omid&rankin 'finite rotation analysis ...'
c          CMANE, vol,93, pp. 353-384;1991; eqns. (81) & (a.17)
c
c  input arguments:
c     nnod   (int,sc,in)            - number of nodes
c     feba   (dp,ar(6*nnod),in)     - fe overbar force vector (local and deformable)
c     deba   (dp,ar(6*nnod),in)     - de overbar displacement vector (local and deformable)
c     nr     (int,sc,in)            - dimension for c matrix
c
c  output arguments:
c     cmtr   (dp,ar(nr,nr),out)     - c matrix used to calculate stress stiffnessing
c
#include "impcom.inc"
c
      integer nnod,i,j,k,k1,nr
      double precision
     & cmtr(nr,nr),feba(6*nnod),deba(6*nnod),wk1(3,3),
     & wk2(3,3),wk3(3,3),v_dot,con,con1,con2,con3,con4,con5,
     & eta,nu,wk4(3,3),wk6(3,3),SMALL
      parameter (SMALL = 1.0d-5)
c
      call v_zero (cmtr(1,1),nr*nr)
      j = 3
      do 10 i = 1,nnod
         con = v_dot (feba(j+1),feba(j+1),3)
         con1 = v_dot (deba(j+1),deba(j+1),3)
         con5 = con1
         con = sqrt (con)
         con1 = sqrt (con1)
         if (con.le.SMALL) goto 10
         call skewA3x3 (feba(j+1),wk1(1,1))
         do 31 k = 1,3
            do 32 k1 = 1,3
               cmtr(k+j,k1+j) = - 0.50d0*wk1(k,k1)
  32        continue
  31     continue
         if (con1.le.SMALL) goto 10
         con2 = 0.50d0*con1
         con3 = con1*con1*sin(con2)
         eta = sin(con2) - con2*cos(con2)
         eta = eta/con3
         nu = con1*(con1+sin(con1)) - 8.0d0*sin(con2)*sin(con2)
         con4 = 4.0d0*con3*con3
         nu = nu/con4
         con = v_dot (feba(j+1),deba(j+1),3)
         call AxBt_C (deba(j+1),feba(j+1),wk1(1,1),3,3,3,3,3,1)
         call AxBt_C (deba(j+1),deba(j+1),wk2(1,1),3,3,3,3,3,1)
         do 13 k = 1,3
            wk2(k,k) = wk2(k,k) - con5
  13     continue
         call AxBt_C (wk2(1,1),wk1(1,1),wk3(1,1),3,3,3,3,3,3)
         do 42 k = 1,3
            do 41 k1 = 1,3
               cmtr(k+j,k1+j) = cmtr(k+j,k1+j) + eta*(wk1(k,k1) -
     &                          2.0d0*wk1(k1,k)) + nu*wk3(k,k1)
               if (k.eq.k1) cmtr(k+j,k1+j) = cmtr(k+j,k1+j) + eta*con
  41        continue
  42     continue
         call skewA3x3 (deba(j+1),wk1(1,1))
         do 22 k = 1,3
            do 21 k1 = 1,3
               wk1(k,k1) = -0.50d0*wk1(k,k1) + eta*deba(k+j)*deba(k1+j)
               if (k.eq.k1)
     &             wk1(k,k1) = wk1(k,k1) + (1.0d0-con5*eta)
  21        continue
  22     continue
         do 14 k = 1,3
            call v_move (cmtr(j+1,k+j),wk4(1,k),3)
  14     continue
         call AxB_C (wk4(1,1),wk1(1,1),wk6(1,1),3,3,3,3,3,3)
         do 15 k = 1,3
            call v_move (wk6(1,k),cmtr(j+1,k+j),3)
  15     continue
  10  j = j + 6
      return
      end
