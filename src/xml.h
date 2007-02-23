 /*
 				xml.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	XML logging.
*
*	Last modify:	23/02/2007
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#ifndef _PSF_H_
#include "psf.h"
#endif

/*----------------------------- Internal constants --------------------------*/
#ifndef XSL_URL
#define	XSL_URL	"."
#endif
/*--------------------------------- typedefs --------------------------------*/
/*------------------------------- functions ---------------------------------*/

extern int	write_xml(char *filename),
		write_xml_header(FILE *file),
		write_xml_meta(FILE *file, char *error),
		write_xmlconfigparam(FILE *file, char *name, char *unit,
				char *ucd, char *format);
extern void	end_xml(void),
		init_xml(out_data_struct *out, int nfield),
		write_xmlerror(char *filename, char *error);
