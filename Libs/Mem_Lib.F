c --- uses either ansys or own library for memory management

#if defined(ANSYS_MEM)
c --- Memory allocation interface to C
c     Why C? Intel Fortran does not support Realloc and LongInt addressing
 
      function MemAlloc (iLen, LenD, PtrName)
c --- primary function: get a block of space from mem manager and return pointer
c
c  input arguments:
c     iLen (int,sc,in)           - length of the block (in data elements)
c     PtrName (chr*16,sc,in)     - 16 character name for the Block
c     key (int,sc,in)            - length of memory element
c                                             
c  output arguments:
c     MemAlloc (PTRFTN,sc,out)  - pointer to this data block
c                                             
#include "impcom.inc"
#include "constants.inc"
c
      integer iLen, LenD
      character*(*) PtrName
      PTRFTN MemAlloc, fAnsMemAlloc

c      integer*8 dLen,i
c      pointer (MemAllocD,debug)
c      integer debug(*)

      MemAlloc = fAnsMemAlloc(iLen,LenD,PtrName)

c      MemAllocD = MemAlloc
c      selectcase (LenD)
c          case (W_LenInt)
c                dLen = iLen*2/W_IntPerLong
c          case (W_LenByte)
c                dLen = iLen/4
c          case (W_LenDouble)
c                dLen = iLen*2
c          case (W_LenLongInt) 
c                dLen = iLen*2
c      endselect
      
c      do i = 1,dLen
c         debug(i) = int(i)
c      enddo
          
      return
      end

      function MemAllocL (iLenL, LenD, PtrName)
c --- primary function: get a block of space from mem manager and return pointer
c
c  input arguments:
c     iLenL (long,sc,in)        - length of the block (in data elements)
c     PtrName (chr*16,sc,in)    - 16 character name for the block
c     key (int,sc,in)           - length of memory element
c                                             
c  output arguments:
c     MemAllocL(PTRFTN,sc,out)  - pointer to this data block
c                                           
#include "impcom.inc"
#include "constants.inc"
c
      integer LenD
      integer*8 iLenL
      character*(*) PtrName
      PTRFTN MemAllocL, fAnsMemAllocL

c      integer*8 dLen,i
c      pointer (MemAllocLD,debug)
c      integer debug(*), irand
      
      MemAllocL = fAnsMemAllocL(iLenL,LenD,PtrName)

c      MemAllocLD = MemAllocL
c      selectcase (LenD)
c          case (W_LenInt)
c                dLen = iLenL*2/W_IntPerLong
c          case (W_LenByte)
c                dLen = iLenL/4
c          case (W_LenDouble)
c                dLen = iLenL*2
c          case (W_LenLongInt) 
c                dLen = iLenL*2
c      endselect
      
c      do i = 1,dLen
c         debug(i) = -100
c      enddo
      
      return
      end

      function MemRealloc (memPtr, iLen, LenD, PtrName)
c --- primary function: reallocate a block of space from mem manager and return pointer
c
c  input arguments:
c     memPtr (PTRFTN,sc,in)      - original pointer
c     iLen (int,sc,in)           - length of the block (in data elements)
c     PtrName (chr*16,sc,in)     - 16 character name for the block
c     key (int,sc,in)            - length of memory element
c                                             
c  output arguments:
c     MemRealloc (PTRFTN,sc,out) - pointer to this data block
c                                             
#include "impcom.inc"
#include "constants.inc"
c
      integer iLen, LenD
      character*(*) PtrName
      PTRFTN MemRealloc, memPtr, fAnsMemRealloc

      MemRealloc = fAnsMemRealloc(memPtr,iLen,LenD,PtrName)

      return
      end

      function MemReallocL (memPtr, iLenL, LenD, PtrName)
c --- primary function: reallocate a block of space from mem manager and return pointer
c
c  input arguments:
c     memPtr (PTRFTN,sc,in)      - original pointer
c     iLen (long,sc,in)          - length of the block (in data elements)
c     PtrName (chr*16,sc,in)     - 16 character name for the Block
c     key (int,sc,in)            - length of memory element
c                                             
c  output arguments:
c     MemRealloc (PTRFTN,sc,out) - pointer to this data block
c                                             
#include "impcom.inc"
#include "constants.inc"
c
      integer*8 iLenL
      integer LenD
      character*(*) PtrName
      PTRFTN MemReallocL, memPtr, fAnsMemReallocL

      MemReallocL = fAnsMemReallocL(memPtr,iLenL,LenD,PtrName)

      return
      end

      subroutine MemFree (memPtr)
c --- primary function: free a data block, given a pointer
c
c  input arguments:
c     ptr   (PTRFTN,sc,inout)  - pointer for this block
c                                             
c  output arguments:  
c     ptr   (PTRFTN,sc,inout)  - pointer will be set to zero
c                                            
#include "impcom.inc"
c
      PTRFTN memPtr

      call fAnsMemFree(memPtr)
      
      return
      end
#else
c --- Memory allocation interface to C
c     Why C? Intel Fortran does not support Realloc and LongInt addressing
 
      function MemAlloc (iLen, LenD, PtrName)
c --- primary function: get a block of space from mem manager and return pointer
c
c  input arguments:
c     iLen (int,sc,in)           - length of the block (in data elements)
c     PtrName (chr*16,sc,in)     - 16 character name for the block
c     LenD (int,sc,in)           - length of memory element
c                                             
c  output arguments:
c     MemAlloc (PTRFTN,sc,out)   - pointer to this data block
c                                             
#include "impcom.inc"
#include "constants.inc"
c
      integer iLen, LenD
      character*(*) PtrName
      PTRFTN MemAlloc, f_malloc, size

      MemAlloc = PTRFTNNULL
      
      if (iLen .gt. 0) then
          size = iLen*LenD
          MemAlloc = f_malloc(size)
c          call MemDebug(MemAlloc,1,PtrName)
      endif

      return
      end

      function MemAllocL (iLenL, LenD, PtrName)
c --- primary function: get a block of space from mem manager and return pointer
c
c  input arguments:
c     iLenL (long,sc,in)        - length of the block (in data elements)
c     PtrName (chr*16,sc,in)    - 16 character name for the block
c     LenD (int,sc,in)          - length of memory element
c                                             
c  output arguments:
c     MemAllocL(PTRFTN,sc,out)  -  pointer to this data block
c                                           
#include "impcom.inc"
#include "constants.inc"
c
      integer LenD
      integer*8 iLenL
      character*(*) PtrName
      PTRFTN MemAllocL, f_malloc, size

      MemAllocL = PTRFTNNULL
      
      if (iLenL .gt. 0) then
          size = iLenL*LenD
          MemAllocL = f_malloc(iLenL*LenD)
c          call MemDebug(MemAllocL,1,PtrName)
      endif

      return
      end

      function MemRealloc (memPtr, iLen, LenD, PtrName)
c --- primary function: reallocate a block of space from mem manager and return pointer
c
c  input arguments:
c     memPtr (PTRFTN,sc,in)      - original pointer
c     iLen (int,sc,in)           - length of the block (in data elements)
c     PtrName (chr*16,sc,in)     - 16 character name for the block
c     LenD (int,sc,in)           - length of memory element
c                                             
c  output arguments:
c     MemRealloc (PTRFTN,sc,out) - pointer to this data block
c                                             
#include "impcom.inc"
#include "constants.inc"
c
      integer iLen, LenD
      character*(*) PtrName
      PTRFTN MemAlloc, MemRealloc, f_realloc, memPtr, size

      MemRealloc = PTRFTNNULL

c --- check for null
      if (memPtr .eq. PTRFTNNULL) then
          MemRealloc=MemAlloc(iLen,LenD,PtrName)
          goto 999
      endif

c --- make resize to 0 act like a MemFree
      if (iLen .eq. 0) then
          call MemFree(memPtr)
          goto 999
      endif
      
      if (iLen .gt. 0) then
          size = iLen*LenD
          MemRealloc = f_realloc(memPtr,size)
      endif

999   return
      end

      function MemReallocL (memPtr, iLenL, LenD, PtrName)
c --- primary function: reallocate a block of space from mem manager and return pointer
c
c  input arguments:
c     memPtr (PTRFTN,sc,in)      - original pointer
c     iLen (long,sc,in)          - length of the block (in data elements)
c     PtrName (chr*16,sc,in)     - 16 character name for the Block
c     LenD (int,sc,in)           - length of memory element
c                                             
c  output arguments:
c     MemRealloc (PTRFTN,sc,out) - pointer to this data block -- needs to be
c                                     tied to a local variable in the calling 
c                                     routine
c                                             
#include "impcom.inc"
#include "constants.inc"
c
      integer*8 iLenL
      integer LenD
      character*(*) PtrName
      PTRFTN MemAllocL, MemReallocL, f_realloc, memPtr, size

      MemReallocL = PTRFTNNULL

c --- check for null
      if (memPtr .eq. PTRFTNNULL) then
          MemReallocL=MemAllocL(iLenL,LenD,PtrName)
          goto 999
      endif

c --- make resize to 0 act like a MemFree
      if (iLenL .eq. 0) then
          call MemFree(memPtr)
          goto 999
      endif
      
      if (iLenL .gt. 0) then
          size = iLenL*LenD
          MemReallocL = f_realloc(memPtr,size)
      endif

999   return
      end

      subroutine MemFree (memPtr)
c --- primary function: free a data block, given a pointer
c
c  input arguments:
c     ptr   (PTRFTN,sc,inout)  - pointer for this block
c                                             
c  output arguments:  
c     ptr   (PTRFTN,sc,inout)  - pointer will be set to zero
c                                            
#include "impcom.inc"
c
      PTRFTN memPtr

c --- freeing a null pointer is harmless
      if (memPtr .eq. PTRFTNNULL) goto 999

      call f_free(memPtr)
c      call MemDebug(MemAlloc,2,'')
      memPtr = PTRFTNNULL

 999  return
      end
      
      subroutine MemDebug (memPtr, aType, PtrName)
c --- primary function: debug
c
#include "impcom.inc"
#include "constants.inc"
c
      PTRFTN memPtr
      integer aType
      character*(*) PtrName
     
      call lockSMP (W_LockText)
      
      write (8,*) memPtr, aType
      if (aType .eq. 1) then
          write (9,*) memPtr, PtrName
      endif
      
      call unlockSMP (W_LockText)

      return
      end
#endif
      