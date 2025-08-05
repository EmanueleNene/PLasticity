      integer           P_Reorder
	  integer*8         P_pt(64)
	  PTRFTN            P_iperm, P_invp
      common /solver_info/
     &             P_Reorder, P_pt, P_iperm, P_invp