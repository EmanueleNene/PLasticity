      subroutine MakeSvrIndex (nelem,fileIndex,SvrLengL)
c --- primary function: define the SVR file indices
c
c  input arguments:
c     fpstartL     (intL,sc,in)    - file pointer at start of data (ie, after headers)
c     nelem        (int,sc,in)     - number of elements
c
c  output arguments:
c     fileIndex    (int,ar(nelem,2),out) - two-word representation of the
c                                          file pointer for each element
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "SVR_file.inc"
c
      PTRFTN    MemAlloc
      integer   elmget, etyget

      integer   nelem, fileIndex(nelem,2)
      integer*8 lcfileL, SvrLengL
      integer   nodes(W_NdMax),elChar(IELCSZ),elData(EL_DIM),
     &          elord,elem,mat,ityp,nnod,i,
     &          SVRSize,fileSize

#if defined(DEBUG)
      call SubBeg('MakeSvrIndex')
#endif

      lcfileL = 0

c --- loop through each element
      elChar(IETYP) = 0
      do elord = 1,nelem
         elem = M_ElemOrd(elord)

         nnod = elmget(elem,elData(1),nodes(1))

         if (nnod.le.0) cycle

c --- get the element characteristics
         mat  = elData(EL_MAT)
         ityp = elData(EL_TYPE)
         if (ityp.ne.elChar(IETYP)) then
             i = etyget(ityp,elChar(1))
             elChar(IETYP) = ityp
         endif

c --- get the size written to the files
         call CalSVRRecLen (elChar(1),elData(1),mat,SVRSize)
         fileSize = SVRSize

c --- set the pointer
         if (fileSize.gt.0) then
c --- this element puts out some data
             call MakeShortInt (lcfileL,
     &                          fileIndex(elord,1),fileIndex(elord,2))
         else
c --- this element puts out no data
             fileIndex(elord,1) = 0
             fileIndex(elord,2) = 0
         endif

c --- update the file pointer
         lcfileL = lcfileL + fileSize
      enddo

      SvrLengL = lcfileL + 1
      
#if defined(DEBUG)
      call SubEnd('MakeSvrIndex')
#endif

      return
      end

