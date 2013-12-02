/*
*				define.h
*
* Global definitions.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 1997-2013 Emmanuel Bertin -- IAP/CNRS/UPMC
*
*	Author:			Emmanuel Bertin (IAP)
*
*	License:		GNU General Public License
*
*	PSFEx is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
* 	(at your option) any later version.
*	PSFEx is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with PSFEx.  If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		02/12/2013
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/* Check if we are using a configure script here */
#ifndef HAVE_CONFIG_H
#define		VERSION		"3.x"
#define		DATE		"2009-03-30"
#endif
/*------------------------ what, who, when and where ------------------------*/

#define         BANNER          "PSFEx"
#define         EXECUTABLE      "psfex"
#define		MYVERSION	VERSION
#define         COPYRIGHT       "2011-2013 IAP/CNRS/UPMC"
#define		DISCLAIMER	BANNER " comes with ABSOLUTELY NO WARRANTY\n" \
		"You may redistribute copies of " BANNER "\n" \
		"under the terms of the GNU General Public License."
#define		AUTHORS		"Emmanuel BERTIN <bertin@iap.fr>"
#define		WEBSITE		"http://astromatic.net/software/psfex"
#define		INSTITUTE	"IAP http://www.iap.fr"

/*----------------------------- Internal constants --------------------------*/
#define		OUTPUT		stderr		/* where all msgs are sent */
#define		BIG		1e+30		/* a huge number */
#define		MAXCHAR		512		/* max. number of characters */
#define		MAXFILE		32768		/* max number of input files */
#define		ALL_EXTENSIONS	(-1)		/* code for all MEF extensions*/

/*------------ Set defines according to machine's specificities -------------*/
#if _LARGEFILE_SOURCE
#define	FSEEKO	fseeko
#define	FTELLO	ftello
#else
#define	FSEEKO	fseek
#define	FTELLO	ftell
#endif

/*---------------------------- return messages ------------------------------*/

#define		RETURN_OK		0
#define		RETURN_ERROR		(-1)
#define		RETURN_FATAL_ERROR	(-2)

/*--------------------- in case of missing constants ------------------------*/

#ifndef PI
#define PI		3.1415926535898	/* never met before? */
#endif

#ifndef         SEEK_SET
#define         SEEK_SET        0
#endif
#ifndef         SEEK_CUR
#define         SEEK_CUR        1
#endif

#ifndef EXIT_SUCCESS
#define 	EXIT_SUCCESS	0
#endif
#ifndef EXIT_FAILURE
#define		EXIT_FAILURE	-1
#endif

/*------------------- a few definitions to read FITS parameters ------------*/

#define	FITSTOF(k, def)	((point = fitsnfind(buf, k, n))? \
					 atof(strncpy(st, &point[10], 70)) \
					:(def))
#define	FITSTOI(k, def)	((point = fitsnfind(buf, k, n))? \
					 atoi(strncpy(st, &point[10], 70)) \
					:(def))

#define	FITSTOS(k, str, def) \
                { if (fitsread(buf,k,str,H_STRING,T_STRING)!= RETURN_OK) \
                    strcpy(str, (def)); \
                }

/*------------------------------- Other Macros -----------------------------*/

#define	DEXP(x)	exp(2.30258509299*(x))	/* 10^x */

#define QFREAD(ptr, size, afile, fname) \
		if (fread(ptr, (size_t)(size), (size_t)1, afile)!=1) \
		  error(EXIT_FAILURE, "*Error* while reading ", fname)

#define QFWRITE(ptr, size, afile, fname) \
		if (fwrite(ptr, (size_t)(size), (size_t)1, afile)!=1) \
		  error(EXIT_FAILURE, "*Error* while writing ", fname)

#define	QFSEEK(afile, offset, pos, fname) \
		if (FSEEKO(afile, (offset), pos)) \
		  error(EXIT_FAILURE,"*Error*: file positioning failed in ", \
			fname)

#define	QFTELL(afile, pos, fname) \
		if ((pos=FTELLO(afile))==-1) \
		  error(EXIT_FAILURE,"*Error*: file position unknown in ", \
			fname)

#define	QFREE(ptr) \
		{free(ptr); \
		ptr = NULL;}

#define	QCALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)calloc((size_t)(nel),sizeof(typ)))) \
		   { \
		   sprintf(gstr, #ptr " (" #nel "=%lld elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		   error(EXIT_FAILURE, "Could not allocate memory for ", gstr);\
                   }; \
                 }

#define	QMALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)malloc((size_t)(nel)*sizeof(typ)))) \
		   { \
		   sprintf(gstr, #ptr " (" #nel "=%lld elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		   error(EXIT_FAILURE, "Could not allocate memory for ", gstr);\
                   }; \
                 }

#define	QREALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)realloc(ptr, (size_t)(nel)*sizeof(typ))))\
		   { \
		   sprintf(gstr, #ptr " (" #nel "=%lld elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		   error(EXIT_FAILURE, "Could not allocate memory for ", gstr);\
                   }; \
                 }

#define QMEMCPY(ptrin, ptrout, typ, nel) \
		{if (ptrin) \
                  {if (!(ptrout = (typ *)malloc((size_t)(nel)*sizeof(typ)))) \
		     { \
		     sprintf(gstr, #ptrout " (" #nel "=%lld elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		     error(EXIT_FAILURE,"Could not allocate memory for ",gstr);\
                     }; \
                   memcpy(ptrout, ptrin, (size_t)(nel)*sizeof(typ)); \
                   }; \
                 }

#define	RINT(x)	(int)(floor(x+0.5))

#define	NPRINTF		if (prefs.verbose_type == NORM) fprintf

#define	NFPRINTF(w,x)	{if (prefs.verbose_type == NORM) \
				fprintf(w, "\33[1M> %s\n\33[1A",x); \
			else if (prefs.verbose_type == FULL) \
				fprintf(w, "%s.\n", x);}

#define	FPRINTF		if (prefs.verbose_type == FULL)	fprintf

#define QPRINTF		if (prefs.verbose_type != QUIET)	fprintf
#define QIPRINTF(w,x)	{if (prefs.verbose_type == NORM) \
				fprintf(w, "\33[7m%s\33[0m\n", x); \
			else if (prefs.verbose_type == LOG) \
				fprintf(w, "%s\n", x);}

#define QBPRINTF(w,x)	{if (prefs.verbose_type == NORM) \
				fprintf(w, "\33[01;31m%s\33[0m\n", x); \
			else if (prefs.verbose_type == LOG) \
			fprintf(w, "%s\n", x);}
