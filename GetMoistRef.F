      function GetMoistRef (matId)
c --- primary function: obtain this material reference moisture
c
c  input arguments:
c        matId         (int,sc,in)       - element material number
c
c  output arguments:
c        GetMoistRef   (dp,func,out)     - the value of the reference moisture
c
#include "impcom.inc"
#include "bfcom.inc"
c
      integer  matId
      double precision GetMoistRef
c

      GetMoistRef = tref

      return
      end
