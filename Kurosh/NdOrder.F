      subroutine NdOrder(Order)
c --- primary function: Gives the order of nodes for hygrostrain transfer
c Input arguments:
c
c output arguments:
c  Order     (dp,arr(3,8),out)
c
#include "impcom.inc"
c

      
      double precision    Order(3,8)
      
      double precision ONE
      parameter (ONE=1.0d0)
      
      Order(1,1) = -ONE
      Order(2,1) = -ONE
      Order(3,1) = -ONE
      
      Order(1,2) =  ONE
      Order(2,2) = -ONE
      Order(3,2) = -ONE
      
      Order(1,3) =  ONE
      Order(2,3) =  ONE
      Order(3,3) = -ONE
      
      Order(1,4) = -ONE
      Order(2,4) =  ONE
      Order(3,4) = -ONE
      
      
      Order(1,5) = -ONE
      Order(2,5) = -ONE
      Order(3,5) =  ONE
      
      Order(1,6) =  ONE
      Order(2,6) = -ONE
      Order(3,6) =  ONE
      
      Order(1,7) =  ONE
      Order(2,7) =  ONE
      Order(3,7) =  ONE
      
      Order(1,8) = -ONE
      Order(2,8) =  ONE
      Order(3,8) =  ONE
      

            
999   return
      end      