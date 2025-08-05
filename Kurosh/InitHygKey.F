      subroutine InitHyg_key()
c --- primary function: Initializes the H_key options
c
c H_key(1)   : Dimension (1D/2D/3D)
c H_key(2)   : Shape-function order (1:linear, 2:quadratic(with mid-node))
c H_key(3-14): element divisions (see figure below)
c H_key(15)  : Method to evaluate strain in 2D (1:use strain at center point of bond, 2: Average strain over bond, 3:Average strain over bond and free segments)
c number of divisions:
c                                            __________________
c                                           /  /           /  /
c                                          /  /           /  /
c                                         /  /           /  /
c                                        /  /           /  /H_key(10)
c                                       /  /   BEAM2   /  /
c                                      /  /           /  /
c                                     /  /           /  /
c                 ___________________/__/___________/__/___________                    _________H_key(4)_________                  _________H_key(3)_________
c                /                  /  /           /  /           /                   |        __H_key(6)_       |                |        _H_key(5)__       |
c               /------------------/--/--H_key(5)-/--/-----------/                    |       |           |      |                |       |           |      |
c              /    BEAM1         /  /      H_key(6)/           H_key(4)              | H_key(13)         | H_key(11)             | H_key(14)         | H_key(12)
c             /------------------/--/-----------/--/-----------/                      |       |___________|      |                |       |___________|      |
c            /__________________/__/___________/__/___________/                       |__________________________|                |__________________________|
c                  H_key(7)    /  /           /  /    H_key(8)                                    BEAM1                                       BEAM2
c                             /  /           /  /
c                            /  /           /  /H_key(9)
c                           /  /           /  /
c                          /__/__H_key(3)_/__/
c


#include "impcom.inc"
#include "beam_hyg.inc"
      
      integer p,n
      

      p = 1        ! number of flaps (0 - none)
      n = 1        ! number of devisions
      H_key(1 ) = 3
      H_key(2 ) = 2
      H_key(3 ) = n
      H_key(4 ) = n
      H_key(5 ) = 0
      H_key(6 ) = 0
      H_key(7 ) = p*n
      H_key(8 ) = p*n
      H_key(9 ) = p*n
      H_key(10) = p*n
      H_key(11) = 1
      H_key(12) = 1
      H_key(13) = 0
      H_key(14) = 0
      H_key(15) = 2
      H_key(16) = 1
      
      H_nGpt = H_key(2) + 1
      
      return
      end