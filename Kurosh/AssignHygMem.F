      subroutine AssignHygMem(pK,pU,pF,pKK,pFF,pKR,pFR,pBB,pWW,
     &                        pCnctvt,pElType)
#include "impcom.inc"
#include "beam_hyg.inc"
      
      integer ProcNum, iProcSMP
      
      integer     c_Dp, c_Int
      integer     pK,pU,pF,pKK,pFF
      integer     pKR,pFR, pBB,pWW
      integer     pCnctvt,pElType
      
      ProcNum = iProcSMP() + 1
      
      c_Dp  = (ProcNum-1)*H_DpSize
      c_Int = (ProcNum-1)*H_IntSize
      
      pK  = c_Dp + 1
      pU  = pK   + H_Dof*H_Dof
      pF  = pU   + H_Dof
      PKK = pF   + H_Dof
      pFF = pKK  + H_key(1)*H_nElNode*H_key(1)*H_nElNode*H_ElemTypK
      pKR = pFF  + H_key(1)*H_nElNode*H_ElemTypF
      pFR = pKR  + (H_Dof-H_nBC)*(H_Dof-H_nBC)
      pBB = pFR  + (H_Dof-H_nBC)
      pWW = pBB  + 18*H_nElNode*H_nGpt*H_nGpt*H_nGpt*H_ElemTypB
      
      pCnctvt  = c_Int + 1
      pElType  = pCnctvt+ H_nElNode*H_nElem
      
      return
      end
      
      
      