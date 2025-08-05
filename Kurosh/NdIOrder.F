      subroutine NdIOrder(Order)
c --- primary function: Gives the order of nodes for hygrostrain transfer
c Input arguments:
c
c output arguments:
c  Order     (int,arr(3,8),out)
c
#include "impcom.inc"
c

      
      integer    Order(3,8)
      
      Order(1,1) = 1
      Order(2,1) = 1
      Order(3,1) = 1
      
      Order(1,2) = 2
      Order(2,2) = 1
      Order(3,2) = 1
      
      Order(1,3) = 2
      Order(2,3) = 2
      Order(3,3) = 1
      
      Order(1,4) = 1
      Order(2,4) = 2
      Order(3,4) = 1
      
      
      Order(1,5) = 1
      Order(2,5) = 1
      Order(3,5) = 2
      
      Order(1,6) = 2
      Order(2,6) = 1
      Order(3,6) = 2
      
      Order(1,7) = 2
      Order(2,7) = 2
      Order(3,7) = 2
      
      Order(1,8) = 1
      Order(2,8) = 2
      Order(3,8) = 2
      

            
999   return
      end      