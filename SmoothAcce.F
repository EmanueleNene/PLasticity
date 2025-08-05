      subroutine SmoothAcce (nr,ls,uelm)
c --- primary function: smooth acceleration and velocity from transient
c
c  input arguments:
c     nr     (int,sc,in)            - number of rows (total dof, this element)
c     ls     (int,ar(nr),in)        - dof list for the element
c
c  output arguments:
c     uelm   (dp,ar(nr),inout)      - solution vector (u, delui, delun, vel, acel)
c                                     with vel or acel smoothed
c
#include "impcom.inc"
#include "solu_info.inc"
#include "dyn_info.inc"
#include "dof_info.inc"
#include "in_mem.inc"
c
      integer  nr,ls(nr),k
      integer  j,inode,idof,NodeExt,kdof,RotDofKey
      double precision uelm(nr,5), alpha,delta,aold,con
      double precision R(3,3), r1(3), r2(3), avec(3), anew, vnew
      equivalence (alpha,I_Consts(2)),(delta,I_Consts(3))
      
#if defined(DEBUG)
      call SubBeg('SmoothAcce')
#endif

c --- adjust needed velocity and acceleration terms dof by dof
      RotDofKey = 4
      do k = 1,nr
c --- 2nd order integration
         if (alpha.gt.1.0d-6) then
c --- implicit acceleration adjustments 
             inode = (ls(k)-1)/D_NumDof + 1 ! node number
             NodeExt = M_BackList(inode)    ! user node number
             inode = (k-1)/D_NumDof + 1     ! node number
             idof = k - (inode-1)*D_NumDof  ! nodal dof
c --- rotational dof
             if (idof.ge.RotDofKey .and. idof.le.RotDofKey+2) then
                 if (idof.eq.RotDofKey) then ! first rotational dof
c --- compute rotation matrix
c --- (assume all rotational dof have int order = 2)
                     call GetRotMat (uelm(k,A_RowIncS),R(1,1))
c --- rotate/transform angular vel and acc
c     from current config. to previous
                     do j = 1,3
                        r1(j) = R(1,j)*uelm(k,A_RowVelo) +
     &                          R(2,j)*uelm(k+1,A_RowVelo) +
     &                          R(3,j)*uelm(k+2,A_RowVelo)
                        r2(j) = R(1,j)*uelm(k,A_RowAcce) +
     &                          R(2,j)*uelm(k+1,A_RowAcce) +
     &                          R(3,j)*uelm(k+2,A_RowAcce)
                     enddo
c --- calculate old accelerations
                     do j = 1,3
                        avec(j) = I_Prm(22)*uelm(k+j-1,A_RowIncS) +
     &                            I_Prm(23)*r1(j) + I_Prm(24)*r2(j)
                     enddo
                 endif

                 kdof = idof - RotDofKey + 1 !rotational dof

c --- now calculate constant average acceleration
c     after transforming old accelerations to current configuration
                 uelm(k,A_RowAcce) = I_Prm(28)*(R(kdof,1)*avec(1) +
     &                                        R(kdof,2)*avec(2) + 
     &                                        R(kdof,3)*avec(3)) +
     &                                        delta*uelm(k,A_RowAcce)
             else          ! translational dof
                 anew = uelm(k,A_RowAcce)
                 vnew = uelm(k,A_RowVelo)
                 aold = I_Prm(22)*uelm(k,A_RowIncS) +
     &                  I_Prm(23)*vnew + I_Prm(24)*anew
c --- now calculate constant average acceleration
                 uelm(k,A_RowAcce) = I_Prm(28)*aold + delta*anew
             endif
         endif
      enddo
#if defined(DEBUG)
      call SubEnd('SmoothAcce')
#endif

      return
      end
