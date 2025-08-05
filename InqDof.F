      subroutine InqDof (idof,node,idir)
c --- primary function: return information about a specified dof
c
#include "impcom.inc"
#include "in_mem.inc"
#include "dof_info.inc"
c
      integer  idof,node,idir

      node = (idof-1)/D_NumDof
      idir = idof - D_NumDof*node
      node = node + 1
      node = M_BackList(node)
      idir = D_DofList(idir)

      return
      end
