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
*	Last modify:	13/02/2009
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/


#ifndef _CONTEXT_H_
#define _CONTEXT_H_

/*--------------------------------- constants -------------------------------*/

#define		MAXCONTEXT		8	/* max. # of context keys */
#define		CONTEXT_KEEPHIDDEN	0
#define		CONTEXT_REMOVEHIDDEN	1

/*--------------------------- structure definitions -------------------------*/

typedef struct context
  {
  char		**name;			/* Context names */
  int		*group;			/* Context groups */
  int		*pcflag;		/* Flags PC contexts */
  int		ncontext;		/* Total number of contexts */
  int		*degree;		/* Group degrees */
  int		ngroup;			/* Number of context groups */
  double	*pc;			/* PC components */
  int		npc;			/* Number of PC components */
  }	contextstruct;

/*-------------------------------- protos -----------------------------------*/

contextstruct	*context_init(char **names, int *group, int ndim, int *degree,
			int ngroup, int pcexflag);


#endif

