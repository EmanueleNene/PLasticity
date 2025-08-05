      subroutine GetDofNumber (node,dir,dof)
c --- primary function: return the dof number for a given node and direction
c
c  input arguments:
c     node     (int,sc,in)         - the node number
c     dir      (int,sc,in)         - the dof direction for this node
c
c  output arguments:
c     dof      (int,sc,out)        - the global DOF number
c
#include "impcom.inc"
#include "in_mem.inc"
#include "dof_info.inc"
c
      integer   vi_find

      integer node,dir,dof, j

c --- find the node dof in the global dof vector
      j = vi_find(D_DofList(1),dir,D_NumDof)
c --- set the global dof number (using the compressed node list)
      dof = (M_FrwdList(node) - 1)*D_NumDof + j

      return
      end
