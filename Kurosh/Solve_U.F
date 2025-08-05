      subroutine Solve_U(K,F,U,KR,FR,BCLoc,FreeLoc)
c --- primary function: solves the FE sustem to find U
#include "impcom.inc"
#include "beam_hyg.inc"
      
      double precision    K(H_Dof,H_Dof)
      double precision    F(H_Dof), U(H_Dof)
      
      integer             m_solve
      
      double precision    FR(H_Dof-H_nBC)

#if defined(PARDISO)
      double precision    KR(H_Dof-H_nBC,H_Dof-H_nBC)
#else
      double precision    KR(H_Dof-H_nBC,H_Dof-H_nBC+1)
#endif
      
      integer             dummy
      integer             BCLoc(H_nBC), FreeLoc(H_Dof-H_nBC)
      
      call v_zero(U(1),H_Dof)
      
#if defined(PARDISO)
      KR = K(FreeLoc,FreeLoc)
      FR = F(FreeLoc)

      call dposv('L',H_Dof-H_nBC,1,KR,H_Dof-H_nBC,
     &           FR,H_Dof-H_nBC,dummy)

      U(FreeLoc) = FR
#else
      KR(:,1:(H_Dof-H_nBC)) = K(FreeLoc,FreeLoc)
      KR(:,1+(H_Dof-H_nBC)) = F(FreeLoc)
      dummy = m_solve(KR,H_Dof-H_nBC,H_Dof-H_nBC,1)
      U(FreeLoc) = KR(:,1+(H_Dof-H_nBC))
#endif      

      return
      end
          
          
      