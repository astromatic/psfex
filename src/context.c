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
*	Last modify:	15/03/2008
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
#include "poly.h"
#include "psf.h"
#include "field.h"

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
VERSION 22/02/2008
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
  QCALLOC(context->pcflag, int, ndim);
  QCALLOC(groupflag, int, ndim);
/* Copy context names and group indices ...*/
/* ... and remove Principal Components if asked to */
  d2=0;
  for (d=0; d<ndim; d++)
    {
    if ((pcflag = !wstrncmp(names[d], "PC?", 80)))
      {
      context->npc++;
      if (!pcexflag)
        context->pcflag[d] = 1;
      }
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
          context->group[d2]=g2+1;
      g2++;
      }
    }

  context->ngroup = g2;
  free(groupflag);

  return context;
  }


/****** context_apply ********************************************************
PROTO	double context_apply(contextstruct *context, psfstruct **psfs, int npsf)
PURPOSE	Find the principal component (the one with the highest eigenvalue) and
	subtract its contribution from the covariance matrix, using the
	iterative "power" method.
INPUT	Covariance matrix,
	output vector,
	Number of principal components.
OUTPUT  Eigenvalue (variance) of the PC.
NOTES   -.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 15/03/2008
 ***/
void context_apply(contextstruct *context, psfstruct *psf,
		fieldstruct **fields, int ext, int ncat)
  {
   psfstruct		*psf2;
   polystruct		*poly, *poly2;
   contextstruct	*context2;
   const double		prime[]={2.0,3.0,5.0,7.0,11.0,13.0,17.0,19.0,23.0,
				29.0,31.0,37.0,41.0,43.0,47.0,53.0,59.0,61.0};
   double		dpos[MAXCONTEXT],
			dval;
   float		*comp, *comp2, *comp2t;
   int			*polycopyflag, *polycopyflagt,
			c,c2, i, n,n2, p, npc, npix,comp2size,ncontext2;

  for (p=0; p<ncat; p++)
    fields[p]->psf[ext] = psf_copy(psf);

/* No PC dependency: the PSF is simply duplicated */
  if (!context->npc)
    return;

  if (sizeof(prime)/sizeof(double) < context->ncontext)
    error(EXIT_FAILURE, "*Internal Error*: ",
		"not enough prime numbers in context_apply()");

/* Prepare component copying process */
/* Use the first PSF as template */
  poly = psf->poly;
  if (poly->ndim != context->ncontext)
    error(EXIT_FAILURE, "*Internal Error*: ",
		"poly->ndim != context->ncontext in context_apply()");

/* Create a new context with the PCs removed */
  context2 = context_init(context->name, context->group, context->ncontext,
			context->degree, context->ngroup, CONTEXT_REMOVEPC);
  ncontext2 = context2->ncontext;
  poly2 = poly_init(context2->group, context2->ncontext, context2->degree,
		context2->ngroup);

/* Identify the polynomial terms that have to be merged using prime numbers */
  c2 = 0;
  for (c=0; c<context->ncontext; c++)
    if (context->pcflag[c])
      dpos[c] = 1.0;
    else
      dpos[c] = prime[c2++];
  poly_func(poly, dpos);
  for (c=0; c<context2->ncontext; c++)
    dpos[c] = prime[c];
  poly_func(poly2, dpos);

/* The polycopyflag matrix flags terms that will be merged */
  QMALLOC(polycopyflag, int, poly2->ncoeff*poly->ncoeff);
  polycopyflagt = polycopyflag;
  for (n2=0; n2<poly2->ncoeff; n2++)
    {
    dval = poly2->basis[n2];
    for (n=0; n<poly->ncoeff; n++)
      *(polycopyflagt++) = (fabs(poly->basis[n] - dval) < 0.1);
    }

/* Merge PSF components for each PSF */
  npix = psf->size[0]*psf->size[1];
  comp2size = npix*poly2->ncoeff;
  npc = context->npc;
  for (p=0; p<ncat; p++)
    {
/*--- Update PC component values */
    c2 = 0;
    for (c=0; c<context->ncontext; c++)
      if (context->pcflag[c])
        dpos[c] = context->pc[p*npc+c2++];
      else
        dpos[c] = 1.0;
    poly_func(poly, dpos);
    QCALLOC(comp2, float, comp2size);
    polycopyflagt = polycopyflag;
    for (n2=0; n2<poly2->ncoeff; n2++)
      for (n=0; n<poly->ncoeff; n++)
        if (*(polycopyflagt++))
          {
          comp = psf->comp + n*npix;
          comp2t = comp2 + n2*npix;
          dval = poly->basis[n];
          for (i=npix; i--;)
            *(comp2t++) += dval**(comp++);
          }
/*-- Replace the new PSF components; 1000000 is just a big number */
    fields[p]->psf[ext] = psf2 = psf_inherit(context2, psf);
    free(psf2->comp);
    psf2->comp = comp2;    
    }

  return;
  }


/****** context_end ***********************************************************
PROTO   void context_end(contextstruct *context)
PURPOSE free memory allocated by a complete context structure.
INPUT   Context structure pointer.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 21/02/2008
*/
void	context_end(contextstruct *context)
  {
   int	d;

  for (d=0; d<context->ncontext; d++)
    free(context->name[d]);
  free(context->name);
  free(context->group);
  free(context->degree);
  free(context->pcflag);
  free(context->pc);

  return;
  }
