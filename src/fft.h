 /*
 				fft.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN, Institut d'Astrophysique de Paris.
*
*	Contents:	Include for psf.c.
*
*	Last modify:	25/08/98
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*---------------------------- Internal constants ---------------------------*/

/*--------------------------- structure definitions -------------------------*/

/*---------------------------------- protos --------------------------------*/
extern void	fourn(float *, unsigned int *, int, int),
		frt2d(float *, float *, int, int, int);

extern float	*fastconv(float *, float *, float *, int, int);

extern int	autoconv(float *, int, int, float *, int, int,
			float *, int, int);

