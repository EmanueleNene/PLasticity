      subroutine DefLsVect (elChar,nodes,ls,nr)
c --- primary function: define the DOF vector for an element (often called)
c
c  input arguments:
c     elChar   (int,ar(*),in)    - the current element type description
c     nodes    (int,ar(*),in)    - the nodes numbers on this element
c
c  output arguments:
c     ls       (int,ar(nr),out)  - the DOF vector for the element
c     nr       (int,sc,out)      - the length of the DOF vector
c
#include "impcom.inc"
#include "in_mem.inc"
#include "constants.inc"
#include "dof_info.inc"
c
      integer  elChar(*),nodes(*),ls(*), j,i,eldof,nnod,
     &         jj,k,kk,np,nr

      eldof = elChar(KDOFS)
      nnod = elChar(NMNDST)
      np = elChar(NMDFPN)
      nr = np*nnod
      call vi_zero (ls(1),nr)
      j = 1
      kk = 0
      do 115 i = 1,7
         if (btest(eldof,i-1)) then
             kk = kk + 1

             if (D_DofList(j) .gt. i) goto 115

 112         if (D_DofList(j) .lt. i) then
                 j = j + 1
                 if (j .gt. D_NumDof) goto 999
                 goto 112
             else
                 jj = kk
                 do k = 1,nnod
                    if (nodes(k).gt.0)
     &                  ls(jj) = (M_FrwdList(nodes(k))-1)*D_NumDof + j
                    jj = jj + np
                 enddo
                 j = j + 1
                 if (j .gt. D_NumDof) goto 999
             endif
         endif
 115  continue

 999  continue


      return
      end
