 /*
 				globals.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	global declarations.
*
*	Last modify:	23/02/2007
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*----------------------- miscellaneous variables ---------------------------*/
/*------------------------------- functions ---------------------------------*/
extern float		fast_median(float *arr, int n),
			hmedian(float *ra, int n);

extern void		makeit(char **incatnames, int ncat),
			write_error(char *msg1, char *msg2);
