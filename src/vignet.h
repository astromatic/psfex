 /*
 				vignet.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include for psf.c.
*
*	Last modify:	06/08/99
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*----------------------------- Internal constants --------------------------*/

#define	INTERPW		6	/* Interpolation function range (x) */
#define	INTERPH		6	/* Interpolation function range (y) */
#define	INTERPLIM	3.0	/* Interpolation limit */
#define	INTERPFAC	3.0	/* Interpolation envelope factor */

#define	INTERPF(x)	(x==0.0?1.0 \
			:(x>INTERPLIM?0.0:(x<-INTERPLIM?0.0 \
			:sin(PI*x)*sin(PI/INTERPFAC*x)/(PI*PI/INTERPFAC*x*x))))

//#define	INTERPF(x)	(fabs(x)>1.0?0.0 : 1 - fabs(x))
//#define	INTERPF(x)	(fabs(x)>0.5? 0.0:1.0)
/*
#define	INTERPF(x)	(x==0.0?1.0 \
			:(x>INTERPLIM?0.0:(x<-INTERPLIM?0.0 \
			:sin(PI*x)*exp(-x*x/4.4)/(PI*x))))
*/
				/* Lanczos approximation */

/*--------------------------------- typedefs --------------------------------*/

typedef  enum {VIGNET_CPY, VIGNET_ADD, VIGNET_SUB, VIGNET_MUL, VIGNET_DIV}
		vigopenum;

/*---------------------------------- protos --------------------------------*/
extern int	vignet_copy(float *pix1, int w1, int h1,
			float *pix2, int w2, int h2, int idx, int idy,
			vigopenum vigop),
		vignet_resample(float *pix1, int w1, int h1, float *pix2,
			int w2, int h2, double dx, double dy, float step2,
			float stepi);
