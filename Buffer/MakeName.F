      subroutine MakeName (fname,fext,fdirin,pname,kerr)
c --- primary function: create and check a file name from a name, extension, and directory
c
#include "impcom.inc"
#include "constants.inc"
#include "file_sys.inc"
c
      integer   Sys_Dir

      character*(B_NameLen) fname
      character*(B_ExtLen)  fext
      character*(B_DirLen)  fdirin, fdir
      character*(B_DirLen)  userdir
      character*(B_PathLen) pname

      integer  namchr(B_PathLen)
      integer  kerr,nchar,i,temchr(B_PathLen),nsave,temusr(B_DirLen),n


      pname = ' '
      kerr = 0
      do i = 1,B_PathLen
         temchr(i) = 32
         namchr(i) = 32
      enddo

c --- DIRECTORY
      fdir = fdirin
      call AsciiUnpack (fdir,temchr(1),B_DirLen)

c --- find last non-blank character
      nchar = 0
      do i = 1,B_DirLen
         namchr(i) = temchr(i)
         if (temchr(i) .gt. W_ChBlank) nchar = i
      enddo

      if (nchar .gt. 0) then

c --- check to see if first character of the directory is a legal character
c     add other checks on the validity of the directory name here also

#if defined (PCWINNT_SYS)
c --- change / to \
         do i = 1,nchar
            if (namchr(i) .eq. W_ChSlash) then
                fdir(i:i) = '\'
                namchr(i) = W_ChBackSlash
            endif
         enddo

c --- and \ at end of directory name if none there
         if (fdir(nchar:nchar) .ne. '\') then
            if (nchar .lt. B_DirLen) then
                nchar = nchar + 1
                fdir(nchar:nchar) = '\'
                namchr(nchar) = W_ChBackSlash
            endif
         endif
#else
c --- change \ to /
         do i = 1,nchar
            if (namchr(i) .eq. W_ChBackSlash) then
                fdir(i:i) = '/'
                namchr(i) = W_ChSlash
            endif
         enddo
c --- and / at end of directory name if none there
         if (fdir(nchar:nchar) .ne. '/') then
            if (nchar .lt. B_DirLen) then
                nchar = nchar + 1
                fdir(nchar:nchar) = '/'
                namchr(nchar) = W_ChSlash
            endif
         endif
#endif

         if (namchr(1).ge.65 .and. namchr(1).ne.W_ChTilda) goto 10
         if (namchr(1).ge.48 .and. namchr(1).le.57) goto 10
         if (namchr(1) .eq. W_ChPeriod) goto 10

c --- if system is unix allow directory to begin with /
#if defined (PCWINNT_SYS)
         if (namchr(1) .eq. W_ChAster) goto 10
#else
         if (namchr(1) .eq. W_ChSlash) goto 10
#endif
         if (namchr(1) .eq. W_ChTilda) then
             kerr = Sys_Dir(userdir)
             if (kerr .eq. 0) then
                 call AsciiUnpack (userdir,temusr(1),B_DirLen)
                 do i = 1,B_DirLen
                    if (temusr(i) .ne. W_ChBlank) n = i
                 enddo
                 if (n + nchar - 2 .gt. B_DirLen) then
                     call erhandler('MakeName',5000,2,
     &               'Directory name is too long.  %/ Directory=
     &                %C.',0.0d0,fdir)
                      kerr = 1
                      goto 999
                 endif
                 do i = 3,nchar
                    namchr(i+n-2) = namchr(i)
                 enddo
                 call vi_move (temusr(1),namchr(1),n)
                 nchar = nchar + n - 2
             endif
             goto 10
         endif
         call erhandler('MakeName',5000,2,
     &                  'Invalid directory.  %/ Directory= %C.'
     &                  ,0.0d0,fdir)
         kerr = 1
         goto 999

      endif


c --- file name
  10  nsave = nchar
      call AsciiUnpack (fname,temchr(1),B_NameLen)
c --- find last non-blank character
      n = 0
      do i = 1,(B_PathLen-nchar)
         namchr(nchar+i) = temchr(i)
         if (temchr(i) .gt. W_ChBlank) n = i
      enddo
      nchar = nchar + n

c --- check to see if first character of the name is a legal character
      if (nchar .le. nsave) then
         if ((namchr(nsave+1) .le. 47) .or.
     &      ((namchr(nsave+1).ge.58).and.(namchr(nsave+1).le.64))) then
            call erhandler ('MakeName',5000,2,
     &      'Invalid file name.  %/ File name= %C.',0.0d0,fname)
            kerr = 1
            goto 999
         endif
      endif


c --- extension
      nchar = nchar + 1
      nsave = nchar
c --- add extension delimeter - will overwrite later if extension blank
      call AsciiUnpack ('.',namchr(nchar),1)
      call AsciiUnpack (fext,temchr(1),8)
c --- find last non-blank character
      n =  0
      do i = 1,8
         namchr(nchar+i) = temchr(i)
         if (temchr(i) .gt. W_ChBlank) n = i
      enddo
      nchar = nchar + n

c --- is extension blank?
      if (nchar .eq. nsave) then
          namchr(nchar) = W_ChBlank
          nchar = nchar - 1
          goto 30
      endif

c --- check to see if first character of the extension is legal
      if (namchr(nsave+1) .lt. W_ChZero) then
          call erhandler('MakeName',5000,2,
     &   'Invalid file extension.  %/ File extension=
     &    %C.',0.0d0,fext)
          kerr = 1
          goto 999
      endif

c --- fill up the name with blanks
  30  do i = nchar+1,B_PathLen
         namchr(i) = W_ChBlank
      enddo

c --- pack into pname
      call AsciiPack (namchr(1),pname,B_PathLen)

  999 continue

      return
      end
