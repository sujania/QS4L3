#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
//#include <sys/mtio.h>
#include <fcntl.h>
#include <errno.h>
//#include <sys/types.h>
#include <unistd.h>


#define F_STR char *
#define F_STRP(x) x
#define LCBUF 100

cgetenv_(name,mess,lmes,lname,lmess)
	int lname, lmess;
	F_STR  name;
	F_STR  mess;
	int *lmes;
{	char buf[LCBUF];
 	char *i, *ptr;
	char *namek, *messk;
        int lnamei , lmessi;
	int k;
	extern char *getenv();
	lnamei=lname;
	lmessi=lmess;
	for ( k=0,namek=F_STRP(name); ((*namek) != ' ') && ((*namek) != '\0') && k < lnamei ; )
	{	
		if(k >= LCBUF)
		{	printf("cgetenv: internal buffer size exceeded:%d\n",k);
			exit(2);
		}
		buf[k]=*namek; k++; namek++;
	}
	if(k >= LCBUF)
	{	printf("cgetenv: internal buffer size exceeded:%d\n",k);
		exit(2);
	}
	buf[k]='\0';
	ptr=getenv(buf);
	if( ptr == (char *)0 )
	{	*lmes=0;
	}
	else

	{
		for ( i=ptr, k=0, messk=F_STRP(mess); (*i) != '\0' ; )
		{	if ( k >= lmessi )
			{	printf("cgetenv: message buffer size exceeded:%d\n",k);
				exit(2);
			}
			(*messk)=(*i); i++; k++; messk++;
		}
		if ( k >= lmessi )
		{	printf("cgetenv: message buffer size exceeded:%d\n",k);
			exit(2);
		}
		(*messk)='\0';
		*lmes = k;
	}
}



copen_(pname, ichan, iopt, ierrno, inew, mode)
	F_STR           pname;
	int           *ichan, *iopt, *ierrno, *inew, *mode;
{
	extern int      errno;
/*	extern int      open(); */
	int             kopt, knew;
	kopt = (int) (*iopt);
	knew = (int) (*inew);
	switch (kopt) {
	case 0:
		kopt = O_RDONLY;
		break;
	case 1:
		kopt = O_WRONLY;
		break;
	case 2:
		kopt = O_RDWR;
		break;
	case 8:
		kopt = O_APPEND;
		break;
	default:
		fprintf(stderr, "copen: unknown option %d", *iopt);
		exit(2);
		break;
	}
	switch (knew) {
	case 0:
		break;
	case 1:
		kopt = kopt | O_EXCL | O_CREAT;
		break;
	case 2:
		kopt = kopt | O_TRUNC | O_CREAT;
		break;
	case 3:
		kopt = kopt | O_TRUNC |  O_EXCL;
		break;
	default:
		fprintf(stderr, "copen: unknown status %d", *inew);
		exit(2);
		break;
	}
	errno = 0;
	*ichan = (int) open(F_STRP(pname), kopt|FNDELAY, (int)(*mode) );
	*ierrno = (int) errno;
}


cclose_(ichan, ires, ierrno)
	int           *ichan, *ires, *ierrno;
{
	extern int      errno;
	extern int      close();
	errno = 0;
	*ires = (int) close((int) (*ichan));
	*ierrno = (int) errno;
}

cread_(ichan, pbuf, pnbyt, ires, ierrno)
	char           *pbuf;
	int           *ichan, *pnbyt, *ires, *ierrno;
{
	extern int      errno;
	errno = 0;
	*ires = (int) read((int) (*ichan), pbuf, (int) (*pnbyt));
	*ierrno = (int) errno;
}

cmkdir_(pname,mode,ires,ierrno)
        int *ires,*ierrno,*mode;
        F_STR pname;
{

  //extern int errno, mkdir();
        errno=0;
        *ires=(int) mkdir(F_STRP(pname),(int)(*mode));
        *ierrno=(int)errno;
}

/* cusleep_(psec) */
/* 	int           *psec; */
/* { */
/* 	int             idum; */
/* 	idum = usleep((unsigned) (*psec)); */
/* } */

/* cmtio_(ichan, iop, icnt, ires, ierrno) */
/* 	int           *ichan, *iop, *icnt, *ires, *ierrno; */
/* { */
/* 	extern int errno, mtio(); */
/* 	errno=0; */
/* 	*ires = mtio(*ichan,*iop,*icnt); */
/* 	*ierrno = errno; */
/* 	if(errno != 0) (void)perror("cmtio"); */
/* } */

/* int mtio(ichan, iop, icnt) */
/* 	int ichan, iop, icnt; */
/* { */
/* 	extern int      errno; */
/* 	int             kopt, rc; */
/* 	struct mtop     magop; */
/* 	struct mtget     magstat; */
/* 	switch (iop) { */
/* 	case 0: */
/* 		kopt = MTWEOF; */
/* 		break; */
/* 	case 1: */
/* 		kopt = MTFSF; */
/* 		break; */
/* 	case 2: */
/* 		kopt = MTBSF; */
/* 		break; */
/* 	case 3: */
/* 		kopt = MTFSR; */
/* 		break; */
/* 	case 4: */
/* 		kopt = MTBSR; */
/* 		break; */
/* 	case 5: */
/* 		kopt = MTREW; */
/* 		break; */
/* 	case 6: */
/* 		kopt = MTOFFL; */
/* 		break; */
/* 	case 7: */
/* 		kopt = MTNOP; */
/* 		break; */
/* 	default: */
/* 		fprintf(stderr, "cmtio: unknown optiion %d", iop); */
/* 		exit(2); */
/* 		break; */
/* 	} */
/* 	errno = 0; */
/* 	magop.mt_op = kopt; */
/* 	magop.mt_count = icnt; */
/* 	rc = ioctl( ichan, MTIOCTOP, (char *)&magop) ; */
/* 	return ( rc ); */
/* } */

void statf(int ichan, int * size, int * istat){
  //Read size (bytes) and mode (read/write) of file ichan
  struct stat buf;
  int             dummy;
  dummy = fstat(ichan, &buf);
  *size =  buf.st_size;
  *istat = buf.st_mode;
}


void cfstat_(int * ichan, int * size, int * istat, int * ierrno){
  // extern int errno;
  (void) statf(*ichan,size,istat);
  *ierrno = errno;
}

clseek_(ichan, offst, iopt, ires, ierrno)
	int           *ichan, *offst, *iopt, *ires, *ierrno;
{
  //extern int      errno /*, lseek()  */;
  //extern long lseek() ;
	int             kopt;
	kopt = (int) (*iopt);
	switch (kopt) {
	case 0:
		kopt = L_SET;
		break;
	case 1:
		kopt = L_INCR;
		break;
	case 2:
		kopt = L_XTND;
		break;
	default:
		fprintf(stderr, "clseek: unknown optiion %d", *iopt);
		exit(2);
		break;
	}
	errno = 0;
	*ires = (long) lseek((int) (*ichan), (*offst), kopt);
	*ierrno = (int) errno;
}



/* coassend_(ichan,command,response,lr,iwt,lencom,lenresp) */
/* int *ichan, *lr, *iwt, lencom, lenresp; */
/* char *command, *response; */
/* { */
/* 	extern int sendoas(); */
/* 	(void)sendoas(*ichan,command,response,lr,*iwt,lencom,lenresp); */
/* } */


/* csleep_(psec) */
/* 	int           *psec; */
/* { */
/* 	int             idum; */
/* 	idum = sleep((unsigned) (*psec)); */
/* } */

cperror_(s)
	F_STR s;
{
	int             idum;
	perror(F_STRP(s));
}

/* sendoas(ichan,command,response,lr,iwt,lencom,lenresp) */
/* int ichan, *lr, iwt, lencom, lenresp; */
/* char *command, *response; */
/* { */
/* } */

/* static int noretry=0; */

/* coasnoretry_() */
/* {	noretry=1; */
/* } */

/* cgethost_(mess,lmes,lmess) */
/* 	int lmess; */
/* 	F_STR mess; */
/* 	int *lmes; */
/* {	extern int gethostname(); */
/* 	int k, lmessi; */
/* 	char *namek; */
/* 	lmessi=lmess; */
/* 	(void) gethostname(F_STRP(mess),lmessi); */
/* 	for ( k=0,namek=F_STRP(mess); ((*namek) != ' ') && ((*namek) != '\0') && (k < 80) ; k++,namek++) */
/* 	{	 */
/* 	} */
/* 	*lmes=k; */

/* } */
     


/* csetoas_(ichan,ierrno) */
/* 	int *ichan, *ierrno; */
/* {	 */
/* 	extern int setupoas(); */
/* 	int jchan; */
/* 	jchan = *ichan; */
/* 	*ierrno = setupoas(jchan); */
/* } */


/* int setupoas(ichan) */
/* 	int ichan; */
/* { */
/* } */
