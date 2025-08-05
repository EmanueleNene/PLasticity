      subroutine BndHygStrn(BeamProps,Ang,ac,bc,dir,ta,tb,
     &                      Xi,Eps,MH,W)
c --- primary function: Finds longitudinal strains at integration points
c                       using 2D or 3D FEM simulation at bond-scale

#include "impcom.inc"
#include "beam_hyg.inc"
      
      double precision    BeamProps(11,2),Ang,Xi(2)
      double precision    ac(3),bc(3),dir(3),ta(3),tb(3)
      
      double precision    Eps(H_Row,H_Col),MH,W(2)
      
      double precision    X(4,2)
      
      double precision    BndEps(8,2)
      
      double precision    CC(6,6,2)
      
      integer     pK,pU,pF,pKK,pFF
      integer     pKR,pFR, pBB, pWW
      integer     pCnctvt,pElType
      integer     i
      
c      double precision    K(H_Dof,H_Dof)
c      double precision    F(H_Dof), U(H_Dof)
c      double precision    KK(H_key(1)*H_nElNode,H_key(1)*H_nElNode,16)
c      double precision    FF(H_key(1)*H_nElNode,16)
c      
c      integer             Connectivity(H_nElNode,H_nElem)
c      integer             ElemType(H_nElem)
c      
 
      call AssignHygMem(pK,pU,pF,pKK,pFF,pKR,pFR,pBB,pWW,
     &                        pCnctvt,pElType)

      call SetupConnectivity(H_IntVars(pCnctvt),H_IntVars(pElType))
      if (H_key(16) .eq. 0) then
          call SetupKKandFF(BeamProps,Ang,H_DpVars(pKK),H_DpVars(pFF))
      elseif (H_key(16) .eq. 1) then
          call SetupBB(BeamProps,Ang,H_DpVars(pBB),H_DpVars(pWW),CC)
          call SetupKK(Ang,H_DpVars(pBB),H_DpVars(pWW),CC,H_DpVars(pKK))
      endif
      
      call Assemble_K_F(BeamProps,ac,bc,dir,ta,Ang,
     &                  H_DpVars(pKK),H_DpVars(pFF),
     &                  H_DpVars(pBB),H_DpVars(pWW),CC,
     &                  H_IntVars(pCnctvt),H_IntVars(pElType),
     &                  H_DpVars(pK),H_DpVars(pF))
      call Solve_U(H_DpVars(pK),H_DpVars(pF),H_DpVars(pU),
     &             H_DpVars(pKR),H_DpVars(pFR),
     &             H_BCData(H_pBCLoc),H_BCData(H_pFreeLoc))
      call FindEffDef(H_DpVars(pU),BeamProps,ac,bc,dir,ta,Ang,
     &                H_IntVars(pCnctvt),BndEps,W)
      call FindHygStrain(BeamProps,Ang,Xi,BndEps,Eps)
c      call FindHygMom2D(H_DpVars(pU),BeamProps,Ang,
c     &                  H_DpVars(pKK),H_DpVars(pFF),
c     &                  H_IntVars(pCnctvt),H_IntVars(pElType),MH)
      
c      call SetupConnectivity(Connectivity,ElemType)
c      call SetupStiffness(BeamProps,Ang,KK,FF)
c      call Assemble_K_F(KK,FF,Connectivity,ElemType, K,F)
c      call Solve_U(K,F,U)
c      call FindEffDef(U,Ang,Connectivity,Def)
c      call FindStrain(BeamProps,Xi,Def,Eps)

      end