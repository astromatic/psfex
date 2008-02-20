/*
 				context.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Type definitions related to contexts
*
*	Last modify:	20/02/2008
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _CONTEXT_H_
#define _CONTEXT_H_

/*--------------------------------- constants -------------------------------*/

#define		CONTEXT_KEEPPC		0
#define		CONTEXT_REMOVEPC	1

/*--------------------------- structure definitions -------------------------*/

typedef struct context
  {
  char		**name;			/* Context names */
  int		*group;			/* Context groups */
  int		ncontext;		/* Total number of dimensions */
  int		*degree;		/* Context degrees */
  int		ngroup;			/* Number of context groups */
  int		pcflag;			/* Set if PC components found */
  }	contextstruct;

/*-------------------------------- protos -----------------------------------*/

contextstruct	*context_init(char **names, int *group, int ndim, int *degree,
	 int ngroup, int pcexflag);

void		context_end(contextstruct *context);

#endif
