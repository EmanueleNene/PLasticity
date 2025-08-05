      integer C_MaxHistNum
      parameter (C_MaxHistNum=5)

      integer 
     & C_Disp,C_Forc,C_nAugm,C_nPen,C_nConA,C_nConF,
     & C_ContChang,C_nDeBo,C_nPenIni,C_ConnUpd,
     & C_MaxIter,C_NSubIni,C_NSubMin,C_NSubMax,
     & C_CutEle,C_CutMat,C_IterCtrl0,C_IterCtrlR,
     & C_IterCtrlG,C_IterCtrlL,C_KeyMono,C_ForcConv,
     & C_ElGap,C_nIterEsti,C_nIterEsti1b,C_nIterEsti2b,
     & C_nFraBo,C_nPivots

      double precision
     & C_SavedInc,C_tGap,C_ForcCnvHist,C_StepFactGap

      common / conv_stat / 
     & C_Disp,C_Forc,C_nAugm,C_nPen,C_nConA,C_nConF,
     & C_ContChang,C_nDeBo,C_nPenIni,C_ConnUpd,
     & C_CutEle,C_CutMat, C_IterCtrl0,C_IterCtrlR,
     & C_IterCtrlG,C_IterCtrlL,C_KeyMono,C_ForcConv(2),
     & C_NSubIni,C_MaxIter(2),C_NSubMin,C_NSubMax,C_SavedInc(5),
     & C_tGap,C_ForcCnvHist(C_MaxHistNum,2),C_StepFactGap,
     & C_ElGap,C_nIterEsti,C_nIterEsti1b,C_nIterEsti2b,C_nFraBo,
     & C_nPivots

c	C_Disp        - status of displacements convergence (not used)
c	C_Forc        - status of out-of-blance forces convergence
c	C_nAugm       - number of augmentations in augmented Lagrangian
c	C_nPen        - # of too many contact penetration
c	C_nConA       - # of contact status change (adhesion)
c	C_nConF       - # of contact status change (friction)
c	C_nDeBo       - # of contact debonding
c	C_nPenIni     - # of big initial penetrations
c	C_ContChang   - key for contact element connectivity change
c	C_ConnUpd     - request update connectivity of contact elements
c	C_NSubIni     - specified (initial) number of substeps to take
c	C_NSubMax     - specified maximum number of substeps to take
c	C_NSubMin     - specified minimum number of substeps to take
c	C_MaxIter     - maximum allowed iterations
c	C_CutEle      - cut key for element distortion 
c	C_CutMat      - cut key for material nonlinearity
c	C_ElGap       - most critical element for gap time predictions
c	C_tGap        - critical gap time/load prediction
c	C_StepFactGap - step factor for gap preditions
c	C_ForcCnvHist - history of force convergence values log (A_ConvCurVal/(A_ConvTol*abs(A_ConvRef)))
c	C_nIterEsti   - estimate for iteration number at convergence
c	C_nIterEsti1b - 1 back estimate for iteration number at convergence
c	C_nIterEsti2b - 2 back estimate of iteration number at convergence
c	C_MaxHistNum  - max # of convergence values used for convergence prediction
c	C_nPivots     - number of negative pivots