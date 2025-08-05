      function Sys_Open (nunit,fname,status_in,access_in,
     &                   form_in,rw_in,posit_in)
c --- primary function: open a file
c
c  input arguments:
c     nunit      (int,sc,in)     - fortran unit number
c     fname      (chr,sc,in)     - character string of file name
c     status_in  (chr,sc,in)     - 'O',old  'N',new  'U',unknown 'R',replace
c     access_in  (chr,sc,in)     - 'D',direct   'S',serial(sequential)
c     form_in    (chr,sc,in)     - 'U',unformatted   'F',formatted
c     rw_in      (chr,sc,in)     - 'R',read only   'W',write only   'B',both
c     posit_in   (chr,sc,in)     - 'B',begining    'E',end
c
c  output arguments:
c     Sys_Open   (int,sc,out)    - 0=normal return,    other=iostat number
c
c
#include "impcom.inc"
#include "constants.inc"
#include "file_sys.inc"
c
      integer  Sys_Inq, GetStringLen

      integer  Sys_Open
      integer  nunit,istat,iacc,iform,irw,ios,namlen


      character*(*)  status_in,access_in,form_in,rw_in,posit_in
      character*1  status,access,form,rw
      character*(*)  fname
      character*12  aform(2)
      character*8  astat(4)
      character*6  posit


      PTRFTN   lposit(2)

      data  aform  / 'UNFORMATTED ','FORMATTED   ' /
      data  astat  / 'OLD     ','NEW     ','UNKNOWN ','REPLACE '/


#if defined(DEBUG)
      call SubBeg('Sys_Open')
#endif

      iform = 2
      irw   = 3
      iacc  = 2
      istat = 3

c --- truncate incoming labels to char*1
      status = status_in
      access = access_in
      form = form_in
      rw = rw_in

      namlen=GetStringLen(fname)

      if (status .eq. 'O') istat = 1
      if (status .eq. 'N') istat = 2
      if (status .eq. 'R') istat = 4
      if (form   .eq. 'U') iform = 1
      if (access .eq. 'D') iacc  = 1
      if (rw     .eq. 'R') irw   = 1
      if (rw     .eq. 'W') irw   = 2
c --- set read status to old
      if (rw     .eq. 'R') istat = 1
c --- convert position to F90 form
      posit = 'REWIND'
      if (posit_in(1:1).eq.'E') posit = 'APPEND'

c --- check if handle already assigned

      if (B_cUnit(nunit).ne.-1) then
          ios = 999
          goto 999
      endif

c --- check for readonly file - do inquire since not using fortran opens
      if (istat .eq. 1) then
          ios = Sys_Inq (fname,'E') - 1
          if (ios .ne. 0) goto 999
      endif

      if (iacc .eq. 2) then                     !.... SEQUENTIAL
          if (iform .eq. 1) then                 !.... UNFORMATTED
              B_PathAndName(nunit) = fname
              call f_namez(B_PathAndName(nunit))
              call f_open (B_PathAndName(nunit),B_cUnit(nunit),
     &                     irw,istat,ios)
              if (ios .ne. 0) goto 999
              lposit(1) = 0
              lposit(2) = 1
              call f_lseek (B_cUnit(nunit),lposit(1),ios)
          else                                    !.... FORMATTED
              if (rw .eq. 'R') then                !.... READ ONLY (PC)
                  open (unit=nunit,file=fname(1:namlen),
     &              status=astat(istat), ACTION='READ',
     &              access='SEQUENTIAL',position=posit,
     &              form=aform(iform),err=999,iostat=ios)
              else                                 !.... READ/WRITE (PC)
                  open (unit=nunit,file=fname(1:namlen),
     &              status=astat(istat),position=posit,
     &              access='SEQUENTIAL',
     &              form=aform(iform),err=999,iostat=ios)
              endif
          endif
      else                                       !.... DIRECT ACCESS
          B_PathAndName(nunit) = fname
          call f_namez(B_PathAndName(nunit))
          call f_open (B_PathAndName(nunit),B_cUnit(nunit),irw,
     &                 istat,ios)
      endif

999   continue
      Sys_Open = ios

#if defined(DEBUG)
      call SubEnd('Sys_Open')
#endif
      return
      end
