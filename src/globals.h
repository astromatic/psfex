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
*	Last modify:	11/03/2008
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*----------------------- miscellaneous variables ---------------------------*/
/*------------------------------- functions ---------------------------------*/
extern float		fast_median(float *arr, int n),
			hmedian(float *ra, int n);

extern void		makeit(void),
			write_error(char *msg1, char *msg2);
