      integer           H_NumBeams, H_Col, H_Row
      parameter        (H_Col=2, H_Row=8)

      integer(8)        H_key(16)
      integer           H_nElem, H_nNode, H_Dof, H_nElNode
      integer           H_ElemTypK,H_ElemTypF,H_ElemTypB

      integer			H_nGpt
      double precision  H_xG(10),H_wG(10)


      integer			H_DpSize,H_IntSize
      double precision  H_DpVars(*)
      integer			H_IntVars(*)
      pointer			(H_pDpVArs,H_DpVars)
      pointer			(H_pIntVArs,H_IntVars)

      integer           H_nBC
      integer			H_BCData(*)
c      double precision    BCVal(3)
      integer           H_pBCLoc,H_pFreeLoc
      pointer			(H_pBCData,H_BCData)

      common    /hyg_vars/ H_NumBeams,H_key,H_nElem,H_nNode,H_Dof, 
     &                     H_nElNode,H_ElemTypK,H_ElemTypF,H_ElemTypB,
     &                     H_nBC, H_pBCLoc, H_pFreeLoc, H_pBCData,
     &                     H_DpSize,H_IntSize, H_pDpVArs, H_pIntVArs,
     &                     H_nGpt, H_xG, H_wG

c    H_NumBeams - number of beam elements
c    H_Col - number of columns in the hygro-strain array