      subroutine CalProjMat(udef,xyz,tr,gam,pte,nnod,nx,nr)
c --- primary function: compute pte the projector matrices for beam element
c
c     compute the pte matrices which will be used to form
c     the geometric stiffness matrix
c     ref: nour-omid & rankin 'finite rotation analysis ...'
c          CMAME, vol. 93, pp. 353-384, 1991, eqns. (33) & (34)
c
c  input arguments:
c     udef   (dp,ar(6*nnod),in)      - total displacement in the global cartesian system
c     xyz    (dp,ar(nx,nnod),in)     - original element coordinates (in global cs) (x,y,z down)
c     tr     (dp,ar(3,3),in)         - ek; rotation matrix transpose
c     nx     (int,sc,in)             - row dimension on xyz
c     gam    (dp,ar(12,3),in)        - gamma matrix from CalGamma routine
c     nnod   (int,sc,in)             - number of nodes
c     nr     (dp,ar(3,3),in)         - total number of variables
c
c  output arguments:
c     pte    (dp,ar(nr,nr),out)      - the element projector matrix (t implies transpose)
c
#include "impcom.inc"
c
      integer nx,nnod,i,j,j1,k,k1,k2,num,num1,num2,nr
      double precision
     & xyz(nx,nnod),udef(nr),urel(3),
     & tr(3,3),tx(3),pte(nr,nr),gam(nr,3),
     & wk1(3,3),wk2(3,3),wk3(3,3),wk4(3,3),wk5(3,3)
c
      num = 6*nnod
      num1 = 3*nr
      num2 = nr*nr
      call v_zero (pte(1,1),num2)
      j = 0
      do 20 i = 1,nnod
         do 11 k = 1,3
            do 12 j1 = 1,num
 12         pte(j1,j+k+3) = -gam(j1,k)
 11      continue
         if (i.ge.2) then
             call v1mv2_v3 (xyz(1,i),xyz(1,1),tx(1),3)
             call v1mv2_v3 (udef(j+1),udef(1),urel(1),3)
             call v1pv2_v1 (tx(1),urel(1),3)
             call Axv1_v1 (tr(1,1),tx(1),3,3)
             call skewA3x3 (tx(1),wk1(1,1))
             j1 = 0
             do 40 k = 1,nnod
                do 50 k1 = 1,3
                   do 60 k2 = 1,3
                      if (k.eq.1) wk4(k1,k2) = -gam(j1+k1+3,k2)
 60                   wk2(k1,k2) = -gam(j1+k1,k2)
 50                continue
                  call AxB_C (wk2(1,1),wk1(1,1),wk3(1,1),3,3,3,3,3,3)
                  if (k.eq.1)
     &                call AxB_C(wk4(1,1),wk1(1,1),wk5(1,1),3,3,3,3,3,3)
                  do 70 k1 = 1,3
                     do 80 k2 = 1,3
                        if (k.eq.1) pte(j1+k1+3,j+k2) = wk5(k1,k2)
 80                     pte(j1+k1,j+k2) = wk3(k1,k2)
 70                 continue
                    j1 = j1 + 6
 40              continue
         endif
         j = j + 6
 20   continue
      do 30 k = 1,num
         pte(k,k) = pte(k,k) + 1.0d0
 30   continue
c
      return
      end
