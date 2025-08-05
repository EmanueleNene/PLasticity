      subroutine BackSubst (nsize,disp)
c --- primary function: calculate displacement solution
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "dof_info.inc"

      integer   cpget, vi_find
      PTRFTN    MemAlloc

      integer  nsize, nodrl1, nodrl,
     &         i, nsm, jj, j, inode, ii, idf, kk

      double precision   disp(nsize)


      pointer (piWork,Work)
      integer Work(*)


#if defined(DEBUG)
      call SubBeg('BackSubst')
#endif

      if (D_NumCp .gt. 0) then
c --- check for coupled degrees of freedom
          piWork = MemAlloc(D_LenCp,W_LenInt,'piWork          ')
          do 210 i = 1,D_NumCp
             nsm = cpget(i,Work(1))
             j = Work(nsm+1)
             idf = vi_find(D_DofList(1),j,D_NumDof)
             if (idf .le. 0) goto 210
c --- search for used master node
             do 100 j = 1,nsm
                nodrl1 = Work(j)
                if (nodrl1 .gt. D_NumNd) cycle
                inode = M_FrwdList(nodrl1)
c --- check for used node
                if (inode .gt. 0) then
                    kk = (inode-1)*D_NumDof + idf
                    if (disp(kk) .ne. HUGE) goto 110
                endif
 100         continue
c --- none of the nodes are used
             goto 210
c --- node j is a valid master
 110         do 200 ii = 1,nsm
c --- skip master node
                if (ii .eq. j) goto 200
                nodrl = Work(ii)
                if (nodrl .gt. D_NumNd) cycle
                inode = M_FrwdList(nodrl)
                if (inode .le. 0) goto 200
                jj = (inode-1)*D_NumDof + idf
                disp(jj) = disp(kk)
 200         continue
 210     continue
         call MemFree (piWork)
      endif


#if defined(DEBUG)
      call SubEnd('BackSubst')
#endif

      return
      end
