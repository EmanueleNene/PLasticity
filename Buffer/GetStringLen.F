      function GetStringLen (chrstr)
c --- primary function: get length of string
c
c  input arguments:
c     chrstr   (chr,sc,in)        - string whose nonblank length is to be found.
c
c  output arguments:
c     GetStringLen   (int,sc,out)  - actual number of characters in the string
c
#include "impcom.inc"
c
      integer  GetStringLen
      integer  i, lenmax
      character*(*)  chrstr

c --- set function to zero - use intrinsic len function to find max string length
      GetStringLen = 0
      lenmax = len(chrstr)

c --- find first nonblank character from end - note it may be zero
      do 10 i = lenmax,1,-1
         if (ichar(chrstr(i:i)).eq.0) goto 10
         if (chrstr(i:i) .eq. ' ') goto 10
         GetStringLen = i
         goto 999
 10   continue
c
 999  continue
c
      return
      end
