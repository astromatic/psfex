  /*
 				pc.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP, Leiden observatory & ESO)
*
*	Contents:	Stuff related to Principal Component Analysis (PCA).
*
*	Last modify:	31/10/2003
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"types.h"
#include	"globals.h"
#include	"fits/fitscat.h"
#include	"prefs.h"
#include	"fft.h"
#include	"sample.h"
#include	"poly.h"
#include	"psf.h"
#include	"vignet.h"

/****** pc_end ***************************************************************
PROTO   void pc_end(pcstruct *pc)
PURPOSE Free a PC structure and everything it contains.
INPUT   pcstruct pointer.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 17/07/99
 ***/
void	pc_end(pcstruct *pc)
  {
   int	i;

  free(pc->comp);
  free(pc->size);
  free(pc->mx2);
  free(pc->my2);
  free(pc->mxy);
  free(pc->flux);
  free(pc->bt);
  if (pc->code)
    {
    free(pc->code->pc);
    for (i=0; i<pc->code->nparam;i++)
      free(pc->code->param[i]);
    free(pc->code->param);
    free(pc->code);
    }
  free(pc);

  return;
  }


/******************************** pc_convolve ********************************/
/*
Convolve the original PCs with the PSF components.
*/
pcstruct	*pc_convolve(pcstruct *pc, psfstruct *psf)
  {
   pcstruct	*pcc;
   codestruct	*code;
   float	*inputpsf,*inputpc,*output;
   int		i,c,n, ncoeff, npc, npix1,npix2,npix3, ncode,nparam;
/* Allocate memory for the convolved PC structure */
  QCALLOC(pcc, pcstruct, 1);
  pcc->dim = 2+1+1;	/* 1 for the PSF components and 1 for the PCs */
  QMALLOC(pcc->size, int, pcc->dim);

  npc = pc->npc;
  ncoeff = psf->poly->ncoeff;
  pcc->npix = pcc->size[0] = psf->size[0];
  npix1 = (pcc->npix *= (pcc->size[1] = psf->size[1]));
  npix2 = (pcc->npix *= (pcc->size[2] = psf->size[2]));
  pcc->npix *= (pcc->size[3] = npc);
  pcc->npc = npc;
  QCALLOC(pcc->comp, float, pcc->npix);
  QMEMCPY(pc->mx2, pcc->mx2, double, npc);
  QMEMCPY(pc->my2, pcc->my2, double, npc);
  QMEMCPY(pc->mxy, pcc->mxy, double, npc);
  QMEMCPY(pc->flux, pcc->flux, double, npc);
  QMEMCPY(pc->bt, pcc->bt, double, npc);
  if (pc->code)
    {
    QMEMCPY(pc->code, pcc->code, codestruct, 1);
    code = pcc->code;
    ncode = code->ncode;
    nparam = code->nparam;
    QMEMCPY(pc->code->pc, code->pc, float, ncode*npc);
    QMALLOC(code->param, float *, nparam);
    for (i=0; i<nparam; i++)
      QMEMCPY(pc->code->param[i], code->param[i], float, ncode);
    }

  inputpsf = psf->comp;
  npix3 = pc->size[0]*pc->size[1];
  for (c=0; c<ncoeff; c++, inputpsf+=npix1)
    {
    output = pcc->comp+c*npix1;
    inputpc = pc->comp;
    for (n=0; n<npc; n++, inputpc+=npix3,output+=npix2)
      autoconv(inputpc, pc->size[0],pc->size[1],
		n?NULL:inputpsf, psf->size[0],psf->size[1],
		output, pcc->size[0],pcc->size[1]);
    }

/* Free memory */
  autoconv(NULL,0,0,NULL,0,0,NULL,0,0);

  return pcc;
  }


/******************************** pc_orthogon ********************************/
/*
Orthogonalize input PCs (using Gram-Schmidt).
*/
pcstruct	*pc_orthogon(pcstruct *pc2, pcstruct *pc, float pixstep)
  {
   pcstruct	*pco;
   codestruct	*code;
   double	*mat,*vec,*vect,
		val, norm;
   float	*pix,*pixt, *opix,*opixt, *bopix,*bopixt, *ubopix,*ubopixt,
		*tpix,*tpixt, *bpix,*ubpix, *minicomp, *cpc, *cpc2,
		*flux,*dscale;
   int		i,p,p2,p3,n, ncoeff, npc, npix,nbpix,nubpix, width,height,
		ncode, nparam;


  code = (codestruct *)NULL;	/* To avoid gcc -Wall warnings */
  mat = vec = (double *)NULL;	/* ditto */

/* Allocate memory for the convolved PC structure */
  QCALLOC(pco, pcstruct, 1);
  *pco = *pc;	/* the new PCs will have the same attributes */
  npc = pco->npc;
  QMEMCPY(pc->size, pco->size, int, pco->dim);
  QMEMCPY(pc->comp, pco->comp, float, pc->npix); /* Copy original pixels */
  QMEMCPY(pc->mx2, pco->mx2, double, npc);
  QMEMCPY(pc->my2, pco->my2, double, npc);
  QMEMCPY(pc->mxy, pco->mxy, double, npc);
  QMEMCPY(pc->flux, pco->flux, double, npc);
  QMEMCPY(pc->bt, pco->bt, double, npc);
  if (pc->code)
    {
    QMEMCPY(pc->code, pco->code, codestruct, 1);
    code = pco->code;
    ncode = code->ncode;
    nparam = code->nparam;
    QMEMCPY(pc->code->pc, code->pc, float, ncode*npc);
    QMALLOC(code->param, float *, nparam);
    for (i=0; i<nparam; i++)
      QMEMCPY(pc->code->param[i], code->param[i], float, ncode);
    }
  else
    ncode = 0;

  ncoeff = pco->dim>3?pco->size[2]:1;
  width = (int)(pco->size[0]*pixstep+0.4999);
  height = (int)(pco->size[1]*pixstep+0.4999);
  npix = width*height;
  nbpix = pco->size[0]*pco->size[1]*ncoeff;
  nubpix = pc2->size[0]*pc2->size[1];

  QMALLOC(minicomp, float, npix*ncoeff*npc);
  QMALLOC(tpix, float, npix*ncoeff);
  if (ncode)
    {
    QCALLOC(mat, double, npc*npc);
    QMALLOC(vec, double, npc);
    }
  opix = minicomp;
  bopix = pco->comp;
  ubopix = pc2->comp;
  for (p=0; p<npc; p++)
    {
/*-- Resample the current PC */
    vignet_resample(bopix, pco->size[0], pco->size[1],
                opix, width, height, 0.5, 0.5, 1/pixstep, 1.0);
/*-- Store the current vector before being overwritten */
    memcpy(tpix, opix, (size_t)(npix*ncoeff)*sizeof(float));
    pix=minicomp;
    bpix=pco->comp;
    ubpix = pc2->comp;
    if (ncode)
      mat[p*npc+p] = 1.0;
    for (p2=0; p2<p; p2++)
      {
      val = 0.0;
      tpixt=tpix;
      pixt=pix;
      for (n=npix; n--;)
        val += *(tpixt++)**(pixt++);
      opixt=opix;
      for (n=npix*ncoeff; n--;)
        *(opixt++) -= val**(pix++);
/*---- Update the full-resolution, convolved PCs */
      bopixt=bopix;
      for (n=nbpix; n--;)
        *(bopixt++) -= val**(bpix++);
/*---- Update the full-resolution, unconvolved PCs */
      ubopixt = ubopix;
      for (n=nubpix; n--;)
        *(ubopixt++) -= val**(ubpix++);
/*---- Update the 2nd order moments */
      pc2->mx2[p] -= val*pc2->mx2[p2];
      pc2->my2[p] -= val*pc2->my2[p2];
      pc2->mxy[p] -= val*pc2->mxy[p2];
/*---- Update the codebook transformation matrix, if available */
      if (ncode)
        for (p3=0;p3<npc;p3++)
          mat[p*npc+p3] -= val*mat[p2*npc+p3];
      }
/*-- Normalize the previous vector */
    norm = 0.0;
    for (opixt=opix, n=npix; n--;)
      {
      val = *(opixt++);
      norm += val*val;
      }
    norm = sqrt(norm);
    for (n=npix*ncoeff; n--;)
      *(opix++) /= (float)norm;
/*-- Update the full-resolution, convolved PCs */
    for (n=nbpix; n--;)
      *(bopix++) /= (float)norm;
/*-- Update the full-resolution, unconvolved PCs */
    for (n=nubpix; n--;)
      *(ubopix++) /= (float)norm;
/*---- Update the 2nd order moments */
    pc2->mx2[p] /= norm;
    pc2->my2[p] /= norm;
    pc2->mxy[p] /= norm;
/*-- Update the codebook transformation matrix, if available */
    if (ncode)
      for (p2=0; p2<npc; p2++)
        mat[p*npc+p2] /= norm;
    }

  bopix = pco->comp;
  ubopix = pc2->comp;
/* Normalize everything with respect to full-resolution convolved PCs */
  for (p=0; p<npc; p++)
    {
    norm = 0.0;
    bopixt = bopix;
    for (n=nbpix; n--;)
      {
      val = *(bopixt++);
      norm += val*val;
      }
    norm = sqrt(norm);
/*-- Update the full-resolution, convolved PCs */
    for (n=nbpix; n--;)
      *(bopix++) /= (float)norm;
/*-- Update the full-resolution, unconvolved PCs */
    for (n=nubpix; n--;)
      *(ubopix++) /= (float)norm;
/*---- Update the 2nd order moments */
    pc2->mx2[p] /= norm;
    pc2->my2[p] /= norm;
    pc2->mxy[p] /= norm;
/*-- Update the codebook transformation matrix, if available */
    if (ncode)
      for (p2=0; p2<npc; p2++)
        mat[p*npc+p2] /= norm;
    }

  if (ncode)
    {
/*-- This may be subject to changes */
    flux = code->param[0];
    dscale = code->param[2];
/*-- Invert the orthogonalization matrix */
    matinv(mat, npc);
/*-- Derive the new flux components */
    cpc = code->pc;
    for (n=ncode; n--;)
      {
/*---- convert the codebook vectors */
      vect = vec;
      norm = 0.0;
      for (p=0; p<npc; p++)
        {
        val = 0.0;
        cpc2 = cpc;
        for (p2=0; p2<npc; p2++)
          val += mat[p+p2*npc]**(cpc2++);
/*---- Compute the norm in the same time */
        *(vect++) = val;
        norm += val*val;
        }
/*---- Normalize while updating the codebook */
      norm = sqrt(norm);
      vect = vec;
      for (p=npc; p--;)
        *(cpc++) = (float)(*(vect++)/norm);
/*---- Update parameters */
      *(flux++) /= (float)norm;
      *(dscale++) *= pixstep;
      }
    }


/* Free memory */
  free(tpix);
  free(minicomp);
  if (ncode)
    {
    free(mat);
    free(vec);
    }

  return pco;
  }


/********************************** pc_load **********************************/
/*
Load the PC data from a FITS file.
*/
pcstruct	*pc_load(char *filename)
  {
   pcstruct	*pc;
   catstruct	*cat;
   tabstruct	*tab;
   keystruct	*key;
   codestruct	*code;
   char		*head, str[80], *ci;
   float	*pc1,*pc2;
   int		i,j, ncode, nparam, npc, pcoffset;

/* Open the cat (well it is not a "cat", but simply a FITS file */
  if (!(cat = read_cat(filename)))
    error(EXIT_FAILURE, "*Error*: PC file not found: ", filename);

/* OK, we now allocate memory for the PC structure itself */
  QCALLOC(pc, pcstruct, 1);

/* Store a short copy of the PC filename */
  if ((ci=strrchr(filename, '/')))
    strcpy(pc->name, ci+1);
  else
    strcpy(pc->name, filename);

  if (!(tab = name_to_tab(cat, "PC_DATA", 0)))
    error(EXIT_FAILURE, "*Error*: PC_DATA table not found in ",
	filename);

/* Load important scalars (which are stored as FITS keywords) */
  head = tab->headbuf;

/* Dimensionality of the PC mask */
  if (fitsread(head, "PCNAXIS", &pc->dim, H_INT, T_LONG) != RETURN_OK)
    goto headerror;
  if (pc->dim<2 || pc->dim>3)
    error(EXIT_FAILURE, "*Error*: wrong dimensionality for the PC "
	"mask in ", filename);
  QMALLOC(pc->size, int, pc->dim);
  for (i=0; i<pc->dim; i++)
    pc->size[i] = 1;
  pc->npix = 1;
  for (i=0; i<pc->dim; i++)
    {
    sprintf(str, "PCAXIS%1d ", i+1);
    if (fitsread(head, str, &pc->size[i], H_INT,T_LONG) != RETURN_OK)
      goto headerror;
    pc->npix *= pc->size[i];
    }

  pc->npc = pc->size[pc->dim-1];

  ncode = 0;
  fitsread(head, "NCODE   ", &ncode, H_INT, T_LONG);
  fitsread(head, "NCODEPAR", &nparam, H_INT, T_LONG);

/* Load the PC mask data */
  key = read_key(tab, "PC_MASK");
  pc->comp = key->ptr;

  key = read_key(tab, "PC_MX2");
  pc->mx2 = key->ptr;

  key = read_key(tab, "PC_MY2");
  pc->my2 = key->ptr;

  key = read_key(tab, "PC_MXY");
  pc->mxy = key->ptr;

  key = read_key(tab, "PC_FLUX");
  pc->flux = key->ptr;

  key = read_key(tab, "PC_BRATIO");
  pc->bt = key->ptr;

  if (ncode)
    {
    QMALLOC(pc->code, codestruct, 1);
    code = pc->code;
    QMALLOC(code->param, float *, nparam);
    QMALLOC(code->parammod, int, nparam);
    code->ncode = ncode;
    code->nparam = nparam;
    key = read_key(tab, "CODE_PC");
    code->pc = (float *)key->ptr;
    for (i=0; i<nparam; i++)
      {
      sprintf(str, "CODE_P%d", i+1);
      key = read_key(tab, str);
      code->param[i] = (float *)key->ptr;
      sprintf(str, "CODE_M%d", i+1);
      fitsread(head, str, &code->parammod[i], H_INT, T_LONG);
      }
    }

/* Prune unwanted PCs */
  if (pc->npc<prefs.pc_npc)
    {
    sprintf(str, "Only %d PC(s) available",pc->npc);
    warning(str, "");
    }
  else if (prefs.pc_npc>0 && pc->npc!=prefs.pc_npc)
    {
    pc->npix /= pc->size[pc->dim-1];
    pcoffset = pc->npc - prefs.pc_npc;
    npc = pc->size[pc->dim-1] = pc->npc = prefs.pc_npc;
    pc->npix *= pc->npc;
    if (ncode)
      {
      code = pc->code;
      pc1 = pc2 = code->pc;
      for (i=ncode; i--; pc2+=pcoffset)
        for (j=npc;j--;)
          *(pc1++) = *(pc2++);
      }
    }

/* But don't touch my arrays!! */
  blank_keys(tab);

  free_cat(&cat, 1);

  return pc;

headerror:
  error(EXIT_FAILURE, "*Error*: Incorrect or obsolete PC data in ", filename);
  return NULL;	/* To avoid gcc -Wall warnings */

  }


/****** matinv ***************************************************************
PROTO   void matinv(double *mat, int nmat)
PURPOSE Invert a matrix using Gauss-Jordan elimination..
INPUT   Pointer to the 2D array,
        Number of matrix elements.
OUTPUT  -.
NOTES   The matrix inversion is done ``in place''..
AUTHOR  E. Bertin (IAP)
VERSION 31/10/2003
 ***/
void	matinv(double *mat, int nmat)
  {
   double	big,dum,pivinv,temp;
   int		*indxc,*indxr,*ipiv,
		i,icol,irow,j,k,l,ll;

#define	SWAP(a,b)	{temp=(a);(a)=(b);(b)=temp;}

  icol = irow = 0; 	/* To avoid gcc -Wall warnings */

/* Allocate memory for the indices */
  QMALLOC(indxc, int, nmat);
  QMALLOC(indxr, int, nmat);
  QMALLOC(ipiv, int, nmat);
  for (j=0;j<nmat;j++)
    ipiv[j]=0;
  for (i=0;i<nmat;i++)
    {
    big=0.0;
    for (j=0;j<nmat;j++)
      if (ipiv[j] != 1)
        for (k=0;k<nmat;k++)
	  {
          if (!ipiv[k])
            {
            if (fabs(mat[j+nmat*k]) >= big)
              {
              big=fabs(mat[j+nmat*k]);
              irow=j;
              icol=k;
              }
            }
          else if (ipiv[k] > 1)
            error(EXIT_FAILURE,"*Error*: Singular Matrix in orthogonalization!"
		, "");
	  }
    ++(ipiv[icol]);
    if (irow != icol)
      for (l=0;l<nmat;l++)
        SWAP(mat[irow+nmat*l],mat[icol+nmat*l]);
    indxr[i]=irow;
    indxc[i]=icol;
    if (mat[icol+nmat*icol] == 0.0)
      error(EXIT_FAILURE, "*Error*: Singular Matrix in orthogonalization!","");
    pivinv = 1.0/mat[icol+nmat*icol];
    mat[icol+nmat*icol] = 1.0;
    for (l=0;l<nmat;l++)
      mat[icol+nmat*l] *= pivinv;
    for (ll=0;ll<nmat;ll++)
      if (ll != icol)
        {
        dum = mat[ll+nmat*icol];
        mat[ll+nmat*icol] = 0.0;
        for (l=0;l<nmat;l++)
          mat[ll+nmat*l] -= mat[icol+nmat*l]*dum;
        }
    }

  for (l=nmat;l--;)
    if (indxr[l] != indxc[l])
      for (k=0;k<nmat;k++)
        SWAP(mat[k+nmat*indxr[l]],mat[k+nmat*indxc[l]]);

  free(ipiv);
  free(indxr);
  free(indxc);

  return;
  }

