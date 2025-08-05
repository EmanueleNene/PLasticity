      subroutine AddGrav (nr,ls,eMass,lvexst,fExt)
c --- primary function:    compute the acceleration load vector for an element
c
c  input arguments:
c     nr       (int,sc,in)         - total number of dof
c     ls       (int,ar(nr),in)     - dof list
c     eMass    (dp,ar(nr,nr),in)   - mass matrix
c
c  output arguments:
c     LVexst   (int,sc,inout)      - key indicating that an element LV exists
c     fExt     (int,sc,inout)      - element load vector, updated for accels
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "dof_info.inc"
#include "acelcm.inc"
c
      integer nr,ls(nr),LVexst,i,node,nodePrev,dof
      double precision
     & eMass(nr,nr),fExt(nr),aclNode(6),acl(W_NrMax)

      call v_zero (aclNode(1),6)
c --- loop thru each DOF
      nodePrev = 0
      do i = 1,nr
c --- see if a new node
         node = (ls(i)-1)/D_NumDof
c --- get the DOF
         dof = D_DofList(ls(i) - D_NumDof*node)
         node = M_BackList(node+1)
         if (node.ne.nodePrev) then
c --- get the acceleration at this node
             call v_move (accel(1),aclNode(1),3)
             nodePrev = node
         endif
c --- move into appropriate DOF slot
         acl(i) = - aclNode(dof)
      enddo

c --- create load vector
      if (LVexst.eq.0) call v_zero (fExt(1),nr)

      call Axv1_v2 (eMass(1,1),acl(1),fExt(1),nr,nr)
      LVexst = 1

      return
      end
