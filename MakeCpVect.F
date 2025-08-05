      subroutine MakeCpVect
c --- primary function: create vector containing coupled DOF
c
#include "impcom.inc"
#include "constants.inc"
#include "dof_info.inc"
#include "in_mem.inc"
c
      integer   cpget, vi_find
      PTRFTN    MemAlloc

      integer   i, nCp, nmax, nTerms, iCp, eDof, iDof, dof, dof1

      pointer (piWorkL,iWork)
      integer  iWork(*)
#if defined(DEBUG)
      call SubBeg('MakeCpVect')
#endif

c --- exit routine if no CP exist
      if (D_NumCp .le. 0) goto 999
c --- allocate memory and initialize coupled DOF list
      M_PtrSlaveCpDofs = MemAlloc(D_ULen,W_LenInt,'M_PtrSlaveCpDofs')
      do i = 1,D_ULen
         M_SlaveCpDofs(i) = i
      enddo

c --- get work space
      nCp = D_NumCp
      nmax = D_LenCp + 1
      piWorkL = MemAlloc(nmax,W_LenInt,'piWorkL         ')

c --- loop over the defined coupling equations
      do iCp = 1,nCp

c --- get the coupling equation information
         nTerms = cpget(iCp,iWork(1))
         if (nTerms .le. 0) cycle
c --- convert from external to internal node and DOF numbering
         do i = 1,nTerms
            iWork(i) = M_FrwdList((iWork(i)))
         enddo
         eDof = iWork(nTerms+1)
         iDof = vi_find(D_DofList(1),eDof,D_NumDof)
         if (iDof .le. 0) cycle

c --- flag coupled DOF
         dof1 = (iWork(1)-1)*D_NumDof + iDof
         M_SlaveCpDofs(dof1) = -iCp
         do i = 2,nTerms
            dof = (iWork(i)-1)*D_NumDof + iDof
            M_SlaveCpDofs(dof) = dof1
         enddo

      enddo

      call MemFree (piWorkL)

  999 continue

#if defined(DEBUG)
      call SubEnd('MakeCpVect')
#endif

      return
      end
