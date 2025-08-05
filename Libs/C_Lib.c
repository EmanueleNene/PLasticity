#include <fcntl.h>
#include <sys/types.h>
#include "computer.h"
#include <string.h>
#if defined (PCWINNT_SYS)                        /* ALL NTs */
#  define UPPERCASENOUNDER_SYS
#  define UPCASE_SYS
#  define SYSTEM_FLOAT_H
#  include <float.h>
#  include <share.h>
#  define Use_memcpy
#endif

typedef long LONG;

#if defined (FULL_64)
typedef long    INT;
typedef unsigned long int  UNSIGNED_INT;
#else
typedef int		INT;
typedef unsigned int  UNSIGNED_INT;
#endif


#if defined(PCWIN64_SYS)
typedef unsigned __int64   PTR;
#elif defined (FULL_64) || defined (SUN64_SYS)
typedef unsigned long int  PTR;
#elif defined (LINUXIA64_SYS) || defined (LINUXOP64_SYS) || defined (LINUXEM64T_SYS)
#if defined(PGI_COMPILER)
typedef unsigned long int  PTR;
#else
typedef unsigned __int64   PTR;
#endif
#else
typedef unsigned int       PTR;
#endif


#include <stdarg.h>
#include <stdlib.h>
#if defined (PCWINNT_SYS)
#include <stdio.h>
#include <errno.h>
#include <sys\stat.h>
#include <io.h>
#else
#include <sys/errno.h>
#endif
#if !defined (PCWINNT_SYS)
#include <sys/uio.h>
#include <sys/file.h>
#endif

#if defined(PCWINNT_SYS)
    /*  C routines called from Fortran  */
#define f_close_   F_CLOSE
#define f_delet_   F_DELET
#define f_namez_   F_NAMEZ
#define f_open_    F_OPEN
#define f_lseek_   F_LSEEK
#define f_read_    F_READ
#define f_write_   F_WRITE
#define f_malloc_  F_MALLOC
#define f_free_	   F_FREE
#define f_realloc_ F_REALLOC
#define v_move_	   V_MOVE
#define vi_move_   VI_MOVE
#endif

#if defined(PTR64)
                                 /* 1234567890123456 */
static const size_t MAX_ADDRESS = 0xFFFFFFFFFFFFFFFF;
#else      
                                 /* 12345678 */
static const size_t MAX_ADDRESS = 0xFFFFFFFF;
#endif

/*
       OPEN A FILE
*/
void f_open_ (namez, handle, irw, istat, err, ilen)

int ilen;
INT *handle, *err, *irw, *istat ;
char namez[];

{

   if (*irw == 1)
  {
                          /*    perm = O_RDONLY;    */
#if defined (PCWINNT_SYS)
   *err = _sopen_s(handle, namez, _O_RDONLY | _O_BINARY, _SH_DENYWR,
                   _S_IREAD);
#else
   *handle = open(namez, O_RDONLY);
#endif
  }
   else
  {
      if (*istat == 2)
     {
                          /*    perm = O_RDWR | O_TRUNC;  */
#if defined (PCWINNT_SYS)
      *err = _sopen_s(handle, namez, _O_RDWR | _O_CREAT | _O_TRUNC | _O_BINARY,
                      _SH_DENYWR, _S_IREAD | _S_IWRITE);
#else
      *handle = open(namez, O_RDWR | O_CREAT | O_TRUNC, 0666);
#endif
     }
      else
                          /*    do not truncate OLD file  */
     {
#if defined (PCWINNT_SYS)
      *err = _sopen_s(handle, namez, _O_RDWR | _O_CREAT | _O_BINARY, _SH_DENYWR,
                      _S_IREAD | _S_IWRITE);
#else
      *handle = open(namez, O_RDWR | O_CREAT, 0666);
#endif
     }
  }

                          /*  handle < 0 is an error                  */
#if !defined (PCWINNT_SYS)
   *err = 0;
#endif
   if (*handle < 0) *err = errno ;
}

/*
      CLOSE A FILE
*/
void f_close_ (handle, err)

INT *handle, *err;

{
#if !defined (PCWINNT_SYS)
   INT close();
#endif

#if defined (PCWINNT_SYS)
   *err = _close (*handle);
#else
   *err = close (*handle);
#endif
       /*  0 is normal  -1 is error  */
   if (*err < 0)
      *err = -errno;
}

/*
      ADD A TRAILING NULL TO A STRING
*/
void f_namez_(namez, nlen)
int nlen;
char namez[];

{
   INT j;

   for (j =nlen-2 ; j >= 0 ; j--)
     {if (namez[j] != ' ')
        {namez[j+1] = '\0';
        break;
        }
     }
}

/*
      POSITION THE FILE SOMEWHERE
*/
void f_lseek_ (handle, loc, err)
#if defined (LINUXIA64_SYS) || defined (LINUXOP64_SYS) || defined (LINUXEM64T_SYS) || defined (SUN64_SYS)
   INT  *handle, loc[2], *err;
#elif defined (PCWINNT_SYS)
   LONG loc[2], *handle;
   INT *err;
#else
   LONG *handle, loc[2], *err;
#endif

{
#if defined (PCWINNT_SYS)
   __int64 location;
   __int64 pos;
#elif defined (LINUXIA32_SYS)
   off_t location;
   off_t pos;
   off_t lseek64();
#else
   LONG lseek();
   LONG pos;
   unsigned long location;
#endif

#if defined (PCWINNT_SYS)
   location = (__int64)loc[0] * (__int64)loc[1];
#elif defined (LINUXIA32_SYS)
   location = (off_t)loc[0] * (off_t)loc[1];
#else
   location = loc[0];
   location = location * loc[1];
#endif
   pos = 0;
#if defined (PCWINNT_SYS)
   pos = _lseeki64(*handle, location, SEEK_SET );
#else
   pos = lseek(*handle, location, 0 );
#endif
   if (pos < 0)
      *err = -errno;
   else
      *err = 0;
}

/*
      WRITE TO THE FILE
*/
void f_write_ (handle, leng, cbuf, err)
INT *handle, *err;
#if !defined (PCWINNT_SYS)
unsigned *leng;
#else
INT *leng;
#endif
char cbuf[];

{
#if defined (PCWINNT_SYS)
   UNSIGNED_INT ileng;
#else
   INT write();
#endif


#if defined (PCWINNT_SYS)
   ileng = (UNSIGNED_INT) *leng;
   *err = _write (*handle, cbuf, ileng);
#else
   *err = write (*handle, cbuf, *leng);
#endif

   if (*err < 0)
      *err = -errno;

}

/*
      READ FROM THE FILE
*/
void f_read_ (handle, leng, cbuf, err)
INT *handle, *err;
#if !defined (PCWINNT_SYS)
unsigned *leng;
#else
INT *leng;
#endif
char cbuf[];

{
#if defined (PCWINNT_SYS)
   UNSIGNED_INT ileng;
#else
   INT read();
#endif


#if defined (PCWINNT_SYS)
   ileng = (UNSIGNED_INT) *leng;
   *err = _read(*handle, cbuf, ileng);
#else
   *err = read(*handle, cbuf, *leng);
#endif

   if (*err < 0)
      *err = -errno;

}
/*
      DELETE THE FILE
*/
void f_delet_ (err,  namez, nlen1)
int nlen1;
char *namez;
INT *err ;

{
#if !defined (PCWINNT_SYS)
   INT unlink();
#endif

/*
 * unlink:
 * Upon successful completion, a value of 0 is returned.
 * Otherwise, a nonzero value is returned and errno is set to
 * indicate the error.
 */
#if !defined (PCWINNT_SYS)
   *err = unlink(namez);
#else
   *err = _unlink(namez);
#endif
   if (*err)
      *err = errno;
}

/*
 *  These are FORTRAN callable versions of malloc, free, and realloc
 */
PTR f_malloc_(PTR *size) {
   return (PTR)malloc((size_t)*size);
}

void f_free_(PTR *ptr) {
   free((void*)*ptr);
}

PTR f_realloc_(PTR *ptr, PTR *size) {
   return (PTR)realloc((void*)*ptr,(size_t)*size);
}

/*
 *  These are FORTRAN callable versions of memspy functions
 */
void v_move_ (PTR *i1, PTR *i2, INT *n)
{
  INT    lth=*n;
  size_t N=lth*sizeof(double);
  memcpy(i2,i1,N);
}
void vi_move_ (PTR *i1, PTR *i2, INT *n)
{
  INT    lth=*n;
  size_t N=lth*sizeof(INT);
  memcpy(i2,i1,N);
}