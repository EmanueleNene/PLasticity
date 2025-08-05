      function user09()
c --- primary function: read in Stefan's files
c
#include "impcom.inc"
#include "constants.inc"
#include "file_sys.inc"
c
      integer user09, Sys_Open, GetWriteUnit, Sys_Close
      integer nxyz, kerr, iott, i, nunit, node, ncount,
     &        elem, ecount, ndpel, nodes(4), valenstal

      double precision xyz(3)

      character*(B_NameLen)  fname
      character*(B_DirLen)   fdir  
      character*(B_PathLen)  pname
      character*32  ch32infun
      parameter (nxyz = 3)

      call Buf_Ini
      
      fname = ch32infun(2)
      call getjobinfo ('CDIR',fdir,i)
      nunit = 9
      ndpel = 2
      iott = GetWriteUnit()
      kerr = 0
     
c --- IMPORT NODES

      call MakeName (fname,'nod     ',fdir,pname,kerr)
      
      kerr = Sys_Open (nunit,pname,'UN','SERIAL','FORM','READ','BEGIN')
      
      if (kerr .ne. 0) goto 998

      ncount = 0
      
  10  read (nunit,*,err=990,end=990)  node,(xyz(i),i=1,nxyz),valenstal
      
      if (node .le. 0) goto 990

      ncount = ncount + 1

c --- store the node      
      call kppxyz (node,0,xyz(1),kerr)

      goto 10

 990  write (iott,2001) ncount
 2001 format (/'Nodes imported:',i8)
 
      kerr = Sys_Close(nunit,'KEEP')

      if (ncount.eq.0) goto 998

c --- IMPORT FIBERS
      call MakeName (fname,'edg     ',fdir,pname,kerr)
     
      kerr = Sys_Open (nunit,pname,'UN','SERIAL','FORM','READ','BEGIN')

      if (kerr .ne. 0) goto 998

      ecount = 0

      
  20  read (nunit,*,err=991,end=991)  elem,(nodes(i),i=1,ndpel)

      if (elem .le. 0) goto 991
      
      if (nodes(1) .eq. nodes(2)) goto 20

  21  ecount = ecount + 1
c --- store the element + 1 zero node
      call lsgbsc (1,nodes(1),nodes(2),1,1.0d0,ecount)
      
      goto 20

 991  write (iott,2002) ecount
 2002 format (/'Segments imported:',i8)

      kerr = Sys_Close(nunit,'KEEP')
            
      user09 = 0
      goto 999

998   continue

      user09 = -999

999   continue

      return
      end