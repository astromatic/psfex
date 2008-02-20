 /*
				context.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Manage observation contexts.
*
*	Last modify:	19/02/2008
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "types.h"
#include "globals.h"
#include "fits/fitscat.h"
#include "prefs.h"
#include "context.h"


/****** context_init *********************************************************
PROTO   contextstruct *context_init(char **names, int *group, int ngroup,
		int *group_deg)
PURPOSE Allocate and initialize a context structure.
INPUT   Pointer to an array of context names,
	Pointer to an array of group indices,
	Number of groups,
	Pointer to an array of dimensions per group,
	Principal Component exclusion flag.
OUTPUT  -.
NOTES   See prefs.h.
AUTHOR  E. Bertin (IAP)
VERSION 19/02/2008
*/
contextstruct	*context_init(char **names, int *group, int ndim, int *degree,
	 int ngroup, int pcexflag)
  {
   contextstruct	*context;
   int			*groupflag,
			d,d2,g,g2, pcflag;

  QCALLOC(context, contextstruct, 1);
  QCALLOC(context->name, char *, ndim);
  QCALLOC(context->group, int, ndim);
  QCALLOC(context->degree, int, ndim);
  QCALLOC(groupflag, int, ndim);
/* Copy context names and group indices ...*/
/* ... and remove Principal Components if asked to */
  d2=0;
  for (d=0; d<ndim; d++)
    {
    pcflag = !wstrncmp(names[d], "PC?", 80);
    context->pcflag |= pcflag;
    if (!pcexflag || !pcflag)
      {
      QMALLOC(context->name[d2], char, 80);
      strncpy(context->name[d2], names[d], 80);
      context->group[d2] = group[d];
      groupflag[group[d]-1]++;
      d2++;
      } 
    }

  context->ncontext = d2;

/* Reorganize groups to remove missing group indices */
  g2 = 0;
  for (g=0; g<ngroup; g++)
    {
    if (groupflag[g])
      {
      context->degree[g2] = degree[g];
      for (d2=0; d2<context->ncontext; d2++)
        if (context->group[d2]==g+1)
          context->group[d2]==g2+1;
      g2++;
      }
    }

  context->ngroup = g2;
  free(groupflag);

  return context;
  }


/****** context_end ***********************************************************
PROTO   void context_end(contextstruct *context)
PURPOSE free memory allocated by a complete context structure.
INPUT   Context structure pointer.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 19/02/2008
*/
void	context_end(contextstruct *context)
  {
   int	d;

  for (d=0; d<context->ncontext; d++)
    free(context->name[d]);

  free(context->group);
  free(context->degree);

  return;
  }
