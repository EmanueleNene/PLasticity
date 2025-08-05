      subroutine DispInGlobal (nr,elChar,nodes,ls,ULen,disp,ndvect,uelm)
c --- primary function: get the displacements for an element in global CS
c
c  input arguments:
c        nr       (int,sc,in)             - number of rows in stiffness matrix
c        elChar   (int,ar(*),in)          - elment type characteristics
c        nodes    (int,ar(*),in)          - the element nodes
c        ls       (int,ar(*),in)          - DOF list for the element
c        ULen     (int,sc,in)             - length of disp array
c        disp     (dp,ar(ULen,ndvect),in) - global disp vectors in nodal cs
c        ndvect   (int,sc,in)             - number of vectors to retrieve
c
c  output arguments:
c        uelm     (dp,ar(nr,ndvect),out) - element disp vectors in global cs
c
#include "impcom.inc"
c
      integer nr,elChar(*),nodes(*),ls(nr),ULen,ndvect

      double precision  disp(ULen,ndvect),uelm(nr,ndvect)
      
c --- get the nodal displacement vectors
      call GetU_Vects (nr,ndvect,ls(1),ULen,disp(1,1),uelm(1,1))
c --- get the pseudovector rotations (Euler parameters) and transfer to global cs
      call GetRot (nr,elChar(1),nodes(1),ndvect,uelm(1,1))


      return
      end
