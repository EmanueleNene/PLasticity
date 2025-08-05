      integer     B_MaxPage,B_MaxBuf
      parameter  (B_MaxPage=8,B_MaxBuf=3)
      integer     B_PathLen,B_ExtLen,B_DirLen,B_NameLen

c --- B_PathLen=B_NameLen(length of filename) + B_ExtLen(length of extension)
c --- + 1(for delimiter) + 1(trailing null) + 2 (pad to 4 byte boundary)
      parameter (B_NameLen=248,B_ExtLen=8,B_DirLen=B_NameLen,
     &           B_PathLen=B_NameLen+B_ExtLen+4)

      integer    B_NumOfPages,   B_UsedPageLen
#if defined(PTR64)
      parameter (B_NumOfPages=B_MaxPage, B_UsedPageLen=65436)
#else
      parameter (B_NumOfPages=B_MaxPage, B_UsedPageLen=65436)
#endif

      integer    B_MaxUnit, B_ByteInWord
      parameter (B_MaxUnit=70, B_ByteInWord=4)
      integer    B_cUnit(B_MaxUnit)
      character*(B_PathLen) B_PathAndName(B_MaxUnit)
      character*(B_NameLen) B_JobName
      character*(B_DirLen) B_Dir      
                                          
      PTRFTN    B_StartPtr
      integer*8 B_StartLoc, B_EndLoc, B_CurPos,
     &          B_MaxLen, B_CurAddr
      integer  B_NumBlock, B_PageLen, B_FortUnit, B_OpenStat,
     &         B_RwStat, B_MaxWrRec, B_CurBlock,
     &         B_CurPage, B_PageNum, B_PageSize,
     &         B_BlockStat, B_StartWord, B_PageAge, B_OffsetStat,
     &         B_PageSizeBit, B_PageFormat, B_AgeCount


      common /file_sys/ 
     & B_StartLoc(B_MaxPage,B_MaxBuf),
     & B_EndLoc(B_MaxPage,B_MaxBuf),
     & B_CurPos(B_MaxBuf), 
     & B_MaxLen(B_MaxBuf),
     & B_CurAddr(B_MaxBuf),
     & B_StartPtr(B_MaxBuf),
     & B_NumBlock(B_MaxBuf),
     & B_PageLen(B_MaxBuf),
     & B_FortUnit(B_MaxBuf),
     & B_OpenStat(B_MaxBuf),
     & B_RwStat(B_MaxBuf),
     & B_MaxWrRec(B_MaxBuf),
     & B_CurBlock(B_MaxBuf),
     & B_CurPage(B_MaxBuf),
     & B_PageNum(B_MaxPage,B_MaxBuf),
     & B_BlockStat(B_MaxPage,B_MaxBuf),
     & B_StartWord(B_MaxPage,B_MaxBuf),
     & B_PageAge(B_MaxPage,B_MaxBuf),
     & B_OffsetStat(B_MaxBuf),
     & B_PageSizeBit(B_MaxBuf),
     & B_PageFormat(B_MaxBuf),
     & B_AgeCount, B_PageSize,
     & B_cUnit, B_PathAndName,
     & B_JobName, B_Dir
     
c	B_MaxPage - maximum number of pages for each file
c	B_MaxBuf - maximum number of file units open

c	B_StartLoc(B_MaxPage,B_MaxBuf) - starting location for each page in memory
c	B_EndLoc(B_MaxPage,B_MaxBuf)   - ending location for each page in memory
c	B_CurPos(B_MaxBuf)- current position on file (32 bit words)
c	B_MaxLen(B_MaxBuf)  - maximum length of file (32 bit words)
c	B_CurAddr(B_MaxBuf)       - base address for current i/o operations
c	B_NumBlock(B_MaxBuf)   - number of blocks for each open file
c	B_PageLen(B_MaxBuf)   - length (32 bit words) for each page
c	B_FortUnit(B_MaxBuf)       - fortran unit number for each file
c	B_OpenStat(B_MaxBuf)       - 0, file not open    1, file open
c	B_RwStat(B_MaxBuf)  - 0,both  1,read  2,write
c	B_MaxWrRec(B_MaxBuf)  - maximum record number written
c	B_CurBlock(B_MaxBuf)   - current block for i/o processing
c	B_CurPage(B_MaxBuf)    - number of in-core pages containing data

c	B_PageNum(B_MaxPage,B_MaxBuf)  - page no. in core page slot
c	B_BlockStat(B_MaxPage,B_MaxBuf)- 0,clean page  1,page written
c	B_StartWord(B_MaxPage,B_MaxBuf)   - starting word for page in core
c	B_PageAge(B_MaxPage,B_MaxBuf)         - age of this page (oldest dumped first)
c	B_OffsetStat(B_MaxBuf)   - 0, relative addressing   1, absolute address
c	B_StartPtr(B_MaxBuf)   - address of first word of buffer
                                                                 
c	B_AgeCount     - counter for setting age
c	B_PageSize    - record size (32 bit integers) for i/o blocks
c	B_PageSizeBit - page size power of 2

c	B_MaxUnit - maximum unit number used
c	B_cUnit     - c file handle for each file (initialize to -1)
c	B_ByteInWord   - number of bytes per integer word (usually 4)
c	B_PathAndName     - the path name of the file

c	B_NameLen     - allowable filename length (excluding extension and directory)
c	B_DirLen     - allowable directory name length
c	B_PathLen     - total allowable length of path name
c	B_NumOfPages   - minimum of 6 for PC is a must!
c	B_UsedPageLen   - page length in integer*4 words (should always be a multiple of 512)
