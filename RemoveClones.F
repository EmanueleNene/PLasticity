      subroutine RemoveClones (isav,irem,nr,matrix)
c --- primary function: remove a duplicate row from element matrix
c
c  input arguments:
c     isav     (int,sc,in)           - row to keep
c     irem     (int,sc,in)           - row to remove
c     nr       (int,sc,in)           - number of rows in matrix
c
c  output arguments:
c     matrix   (int,ar(nr,nr),inout) - modified matrix
c
#include "impcom.inc"
c
      integer  isav,irem,nr
      double precision  matrix(nr,nr)

      call v1pv2_v1i (matrix(isav,1),nr,matrix(irem,1),nr,nr)
      call v_zeroi (matrix(irem,1),nr,nr)
      call v1pv2_v1 (matrix(1,isav),matrix(1,irem),nr)
      call v_zero (matrix(1,irem),nr)

      return
      end
