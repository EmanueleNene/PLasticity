      subroutine PresLoadStiff (nUsrDof,numNodesFace,nDofNode,
     &                          nTranDofNode,unsym,prndv,shpe,
     &                          work,gStiff)
c --- primary function: calculate pressure load stiffness at an integration point
c
c  input arguments:
c     nUsrDof      (int,sc,in)       - number of rows/columns in gStiff matrix
c     numNodesFace (int,sc,in)       - number of nodes on element face (2 or 3)
c     nDofNode     (int,sc,in)       - number of dof per node
c     nTranDofNode (int,sc,in)       - number of translational dof per node
c     unsym        (int,sc,in)       - 0 = symmetric matrix,
c                                      1 = unsymmetric matrix
c     prndv        (dp,ar(3),in)     - pressure*(dL)*normal, where
c                                      dL = geometric weighting factor
c                                      normal = cross product of applied load and edge direction
c     shpe         (dp,ar(3,2),in)   - each row refers to node number,
c                                      first column refers to shape function value,
c                                      second column refers to actual spatial derivative
c     work         (dp,ar(numNodesFace*nDofNode,numNodesFace*nDofNode),in) - scratch space
c
c  output arguments:
c     gStiff       (dp,ar(nUsrDof,nUsrDof),inout)  - pressure load stiffness matrix
c
#include "impcom.inc"
c
      integer nUsrDof,numNodesFace,nDofNode,nTranDofNode,unsym,
     &        ii,jj,kk,ll,i,j,nmat
      double precision 
     &                prndv(3),shpe(3,2),
     &                work(numNodesFace*nDofNode,numNodesFace*nDofNode),
     &                gStiff(nUsrDof,nUsrDof),tmp(3,3),t1
c


#if defined(DEBUG)
      call SubBeg('PresLoadStiff')
#endif

      tmp(1,1) = 0.0d0
      tmp(2,2) = 0.0d0
      tmp(3,3) = 0.0d0
      tmp(1,2) = prndv(3)
      tmp(1,3) =-prndv(2)
      tmp(2,1) =-prndv(3)
      tmp(2,3) = prndv(1)
      tmp(3,1) = prndv(2)
      tmp(3,2) =-prndv(1)
      nmat = numNodesFace*nDofNode
      do i = 1,nmat
         do j = 1,nmat
            work(i,j) = 0.0d0
         enddo
      enddo
      do ii = 1,numNodesFace
         kk = (ii-1)*nDofNode
         do jj = 1,numNodesFace
            ll = (jj-1)*nDofNode
            t1 = shpe(jj,1)*shpe(ii,2)
            do i = 1,nTranDofNode
               do j = 1,nTranDofNode
                  work(kk+i,ll+j)=t1*tmp(i,j)
               enddo
            enddo
         enddo
      enddo
      if (unsym.eq.0) then
          do j=1,nUsrDof
             do i=1,nUsrDof
                gStiff(i,j)=gStiff(i,j)+0.5d0*(work(i,j)+work(j,i))
             enddo
          enddo
      else
          do j=1,nUsrDof
             do i=1,nUsrDof
                gStiff(i,j)=gStiff(i,j)+work(i,j)
             enddo
          enddo
      endif

#if defined(DEBUG)
      call SubEnd('PresLoadStiff')
#endif

      return
      end
