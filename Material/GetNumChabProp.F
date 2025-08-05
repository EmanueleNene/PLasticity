      function GetNumChabProp (keycpL)
c *** primary function: get number of defined models in for Chaboche model
c
c  input arguments:
c     keycpL   (int*8,sc,in)       - key indicating the material and property number
c
c  output arguments:
c     GetChabProp    (int,sc,out)  - NPTS defined by TB command
c
#include "impcom.inc"
c

      integer*8        keycpL
      integer          GetNumChabProp 

      double precision Header(*)
      pointer         (dPtrHeader,Header)
      
c --- get the pointer to the start of the table of property data
      call nlgetPtrFromKey(keycpL,dPtrHeader)
      GetNumChabProp = nint(Header(15))
     
      return
      end
