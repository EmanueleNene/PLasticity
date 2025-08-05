      subroutine GetRot (nr,elChar,nodes,ndvect,uelm)
c --- primary function: get total rotation
c
c  input arguments:
c        nr       (int,sc,in)       - number of rows in stiffness matrix
c        elChar   (int,ar(*),in)    - elment type characteristics
c        nodes    (int,ar(*),in)    - the element nodes
c        ls       (int,ar(*),in)    - DOF list for the element
c        ndvect   (int,sc,in)       - number of vectors to retrieve
c        elord    (int,sc,in)       - current element position
c
c  output arguments:
c        uelm     (dp,ar(nr,ndvect),out) - element disp vectors
c
#include "impcom.inc"
#include "in_mem.inc"
#include "constants.inc"
#include "solu_info.inc"
#include "dyn_info.inc"
c
      integer nr,elChar(*),nodes(*),ndvect,i,j,k,rotloc,nodloc
      double precision  uelm(nr,ndvect),r1(3),r2(3),r3(3),r4(3),R(3,3)

c --- load the pseudovector rotations (Euler parameters) into uelm
      if (iand(elChar(KDOFS),ROTZ).eq.0) goto 999
c --- start of ROTX
      rotloc = 4
c --- put the rotations in storage
      do i = 1,elChar(NMNDST)
         if (nodes(i).gt.0) then
             nodloc = M_FrwdList(nodes(i))
             k = rotloc - 1
             do j = 1,3
                uelm(k+j,A_RowDisp) = M_ERot(j,nodloc)
                if (uelm(k+j,A_RowDisp).eq.HUGE)
     &              uelm(k+j,A_RowDisp) = 0.0d0
                if (uelm(k+j,A_RowDisp).eq.TINY)
     &              uelm(k+j,A_RowDisp) = 0.0d0
             enddo
             rotloc = rotloc + elChar(NMDFPN)
         endif
      enddo

c --- rotate V, A
      if (A_Trans.gt.0 .and. A_KeyFirstPass.eq.0) then
c --- find the start of ROTX
          rotloc = 4
c --- rotate
          do i = 1,elChar(NMNDST)
             if (nodes(i).gt.0) then
c --- process nodes
                 nodloc = M_FrwdList(nodes(i))
                 k = rotloc - 1
c --- rotate V and A vectors from tn-1 by rot this step
                 call GetRotMat (uelm(k+1,A_RowIncS),R(1,1))
c --- 1st recover previous V and A
                 do j = 1,3
                    r1(j) = -I_Prm(33)*uelm(k+j,A_RowVelo) -
     &                       I_Prm(32)*uelm(k+j,A_RowAcce)
                    r2(j) = -I_Prm(20)*uelm(k+j,A_RowVelo) -
     &                       I_Prm(21)*uelm(k+j,A_RowAcce)
                 enddo
c --- now rotate
                 do j = 1,3
                    r3(j) = I_Prm(31)*uelm(k+j,A_RowIncS) +
     &                      R(j,1)*r1(1) + R(j,2)*r1(2) + R(j,3)*r1(3)
                    r4(j) = I_Prm(19)*uelm(k+j,A_RowIncS) +
     &                      R(j,1)*r2(1) + R(j,2)*r2(2) + R(j,3)*r2(3)
                 enddo
                 call v_move (r3(1),uelm(k+1,A_RowVelo),3)
                 call v_move (r4(1),uelm(k+1,A_RowAcce),3)
             endif
             rotloc = rotloc + elChar(NMDFPN)
          enddo
      endif

  999 continue

      return
      end
