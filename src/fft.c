/*
                                  fft.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*       Part of:        PSFEx
*
*       Author:         E.BERTIN (IAP)
*
*       Contents:       Routines dealing with FFT.
*
*       Last modify:    31/10/2003
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
#include "fft.h"

#define	SWAP(a,b)	tempr=(a);(a)=(b);(b)=tempr

/********************************** fourn ************************************/
/*
ndim-dimensional Fast Fourier Transform (from Numerical Recipes in C, p.468).
The 1-dimensional 'data' array is replaced by its Fourier transform:
storage in 'data': (Re[0][0],Im[0][0]); (Re[0][1],Im[0][1]);...
The 1-dimensional 'nn' array contains the sizes in each dimension.
isign = +1 for direct FFT, or -1 for inverse FFT.
Note: float version.
*/
void fourn(float *data, unsigned int *nn, int ndim, int isign)
  {
   int		idim;
   unsigned int	i1,i2,i3,i2rev,ip1,ip2,ip3,ifp1,ifp2;
   unsigned int	ibit,n,nprev,nrem,ntot;
   float	tempi,tempr, *d3, *d3rev, *dk1,*dk2;
   double	theta,wi,wpi,wpr,wr,wtemp;

/* a small trick to begin arrays with nb 1 (a reminiscence of FORTRAN!) */
  data--;

  ntot=1;
  for (idim=0;idim<ndim;idim++)
    ntot *= nn[idim];
  nprev=1;
  for (idim=ndim;idim--;)
    {
    n=nn[idim];
    nrem=ntot/(n*nprev);
    ip1=nprev << 1;
    ip2=ip1*n;
    ip3=ip2*nrem;
    i2rev=1;
    for (i2=1;i2<=ip2;i2+=ip1)
      {
      if (i2 < i2rev)
        {
        for (i1=i2;i1<=i2+ip1-2;i1+=2)
          {
          for (i3=i1;i3<=ip3;i3+=ip2)
            {
            d3rev=i2rev+(d3=data+i3)-i2;
            SWAP(*d3,*d3rev);
            SWAP(*(d3+1),*(d3rev+1));
            }
          }
        }
      ibit=ip2 >> 1;
      while (ibit >= ip1 && i2rev > ibit)
        {
        i2rev -= ibit;
        ibit >>= 1;
        }
      i2rev += ibit;
      }
    ifp1=ip1;
    while (ifp1 < ip2)
      {
      ifp2=ifp1 << 1;
      theta=isign*6.28318530717959/(ifp2/ip1);
      wtemp=sin(0.5*theta);
      wpr = -2.0*wtemp*wtemp;
      wpi=sin(theta);
      wr=1.0;
      wi=0.0;
      for (i3=1;i3<=ifp1;i3+=ip1)
        {
        for (i1=i3;i1<=i3+ip1-2;i1+=2)
          {
          for (i2=i1;i2<=ip3;i2+=ifp2)
            {
            dk2=(dk1=data+i2)+ifp1;
            tempr=wr**dk2 - wi**(dk2+1);
            tempi=wr**(dk2+1) + wi**dk2;
            *(dk2++)=*dk1-tempr;
            *dk2=*(++dk1)-tempi;
            *(dk1--) += tempi;
            *dk1 += tempr;
            }
          }
        wr+=(wtemp=wr)*wpr-wi*wpi;
        wi+=wi*wpr+wtemp*wpi;
        }
      ifp1=ifp2;
      }
    nprev *= n;
    }

  return;
  }


/********************************* autoconv **********************************/
/*
Optimized 2-dimensional FFT convolution, dedicated to 2D real data (images). 
Handles automatically padding to powers of 2 and does not modify input
arrays.
*/
int	autoconv(float *data, int iwidth, int iheight,
		float *mask, int width, int height,
		float *output, int owidth, int oheight)
  {
   static float	*maskbuf, *databuf, *speqmask;
   static int	mwidth, mheight;
   float	*buf1,*buf2;
   int		wmax,hmax, xoffset, hxoffset,hyoffset, hwidth, hheight,
		i,j;

/* Free memory and exit if a NULL output pointer is given */
  if (!output)
    {
    free(maskbuf);
    maskbuf = NULL;
    free(speqmask);
    speqmask = NULL;
    free(databuf);
    databuf = NULL;
    return RETURN_OK;
    }


/* Padding of the mask if a non-zero mask is provided */
  if (mask)
    {
/*-- Search for the biggest image dimensions */
    wmax = width;
    hmax = height;
    if (wmax<iwidth)
      wmax = iwidth;
    if (wmax<owidth)
      wmax = owidth;
    if (hmax<iheight)
      hmax = iheight;
    if (hmax<oheight)
      hmax = oheight;
    mwidth = 2;
    for (i=wmax-1; i>>=1; mwidth<<=1);
    mheight = 2;
    for (i=hmax-1; i>>=1; mheight<<=1);
    free (maskbuf);
    free(speqmask);
    free(databuf);
    speqmask = NULL;
    QCALLOC(maskbuf, float, mwidth*mheight);
/*-- Descrambling of the mask (PSF center put at 1,1) */
    hwidth = width/2;
    hheight = height/2;
    xoffset = mwidth - width;
    buf1 = mask+hwidth+width*hheight;
    buf2 = maskbuf;
    for (j=height-hheight; j--;)
      {
      for (i=width-hwidth; i--;)
        *(buf2++) = *(buf1++);
      buf1 -= width;
      buf2 += xoffset;
      for (i=hwidth; i--;)
        *(buf2++) = *(buf1++);
      buf1 += width;
      }
    buf1 = mask+hwidth;
    buf2 += mwidth*(mheight-height);
    for (j=hheight; j--;)
      {
      for (i=width-hwidth; i--;)
        *(buf2++) = *(buf1++);
      buf1 -= width;
      buf2 += xoffset;
      for (i=hwidth; i--;)
        *(buf2++) = *(buf1++);
      buf1 += width;
      }
    QMALLOC(databuf, float, mwidth*mheight);
    }

/* Padding of the data */
  memset(databuf, 0, (size_t)(mwidth*mheight)*sizeof(float));
  xoffset = mwidth - iwidth;
  hxoffset = (xoffset+1)/2;
  hyoffset = (mheight-iheight+1)/2;
  buf1 = data;
  buf2 = databuf+hyoffset*mwidth+hxoffset;
  for (j=iheight; j--; buf2 += xoffset)
    for (i=iwidth; i--;)
      *(buf2++) = *(buf1++);

  speqmask = fastconv(databuf, maskbuf, speqmask, mwidth, mheight);

/* ``Unpadding'' of the data */
  xoffset = mwidth - owidth;
  hxoffset = (xoffset+1)/2;
  hyoffset = (mheight-oheight+1)/2;
  buf1 = databuf+hyoffset*mwidth+hxoffset;
  buf2 = output;
  for (j=oheight; j--; buf1 += xoffset)
    for (i=owidth; i--;)
      *(buf2++) = *(buf1++);


  return RETURN_OK;
  }

/********************************* fastconv **********************************/
/*
Optimized 2-dimensional FFT convolution, dedicated to 2D real data (images). 
Returns the Nyquist vertical frequency vector (for which memory is allocated).
*/
float *fastconv(float *image, float *mask, float *speqmask,
	int width, int height)
  {
   int		i, nbfreq;
   float	*speqimage, *imap, *maskp,
		real, imag, fac;

/* Allocate temporary memory for storing Nyquist unfolded freq. components */
  QMALLOC(speqimage, float, 2*height);
  if (!speqmask)
    {
    QMALLOC(speqmask, float, 2*height);
/*-- Forward FFT for mask */
    frt2d(mask, speqmask, width, height, 1);
    }
/* Forward FFT for image */
  frt2d(image, speqimage, width, height, 1);
/* Actual convolution (Fourier product) */
  nbfreq = width*height/2.0;
  fac = 1.0/nbfreq;
  imap = image;
  maskp = mask;
  for (i=0; i<nbfreq; i++, maskp+=2)
    {
    real = *imap**maskp-*(imap+1)**(maskp+1);
    imag = *imap**(maskp+1)+*(imap+1)**maskp;
    *(imap++) = fac*real;
    *(imap++) = fac*imag;
    }

  nbfreq = height;
  imap = speqimage;
  maskp = speqmask;
  for (i=0; i<nbfreq; i++, maskp+=2)
    {
    real = *imap**maskp-*(imap+1)**(maskp+1);
    imag = *imap**(maskp+1)+*(imap+1)**maskp;
    *(imap++) = fac*real;
    *(imap++) = fac*imag;
    }

/* Reverse FFT */
  frt2d(image, speqimage, width, height, -1);
  free(speqimage);

  return speqmask;
  }

/*********************************** frt2d ***********************************/
/*
Fast real 2D Fourier transform, based on FFT (inspired from Numerical Recipes
in C, 2nd ed., p.528). 
*/
void frt2d(float *data, float *speq, int width, int height, int isign)
  {
   unsigned int	i2,i3,ii3,j2,j3, nn[2], uwidth, uheight;
   float	*datai2,*datai21, *dataj2,*dataj21, *speqj2,*speqj21,
		c2,h1r,h1i,h2r,h2i;
   double	theta, wi,wpi,wr,wpr,wtemp;

  uwidth = width;
  uheight = height;
  theta = isign*6.28318530717959/height;
  wtemp = sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi = sin(theta);
  nn[0] = uheight;
  nn[1] = uwidth>>1;
  c2 = -0.5*isign;
  if (isign==1)
    {
    fourn(data, nn, 2, isign);
    datai2 = data;
    speqj2 = speq;
    for (i2=0; i2<uheight; i2++, datai2 += uwidth)
      {
      *(speqj2++) = *(datai2++);
      *(speqj2++) = *(datai2--);
      }
    }
  wr = 1.0;
  wi = 0.0;
  for (i3=0,ii3=0; i3<=(uwidth>>2); i3++,ii3+=2)
    {
    if (!i3)
      {
      datai2 = data;
      for (i2=0; i2<uheight; i2++, datai2+=uwidth)
        {
        datai21 = datai2 + 1;
        j2 = i2? ((uheight-i2-1)<<1)+2 : 0;
        speqj2 = speq+j2;
        speqj21 = speqj2+1;
        h1r = 0.5*(*datai2+*speqj2);
        h1i = 0.5*(*datai21-*speqj21);
        h2i = c2*(*datai2-*speqj2);
        h2r = -c2*(*datai21+*speqj21);
        *datai2 = h1r+h2r;
        *datai21 = h1i+h2i;
        *speqj2 = h1r-h2r;
        *speqj21 = h2i-h1i;
        }
      }
    else
      {
      datai2 = data+ii3;
      for (i2=0; i2<uheight; i2++, datai2+=uwidth)
        {
        datai21 = datai2 + 1;
        j2 = i2? uheight-i2 : 0;
        j3 = uwidth+2-((i3+1)<<1);
        dataj2 = data + j2*uwidth + j3;
        dataj21 = dataj2 + 1;
        h1r = 0.5*(*datai2+*dataj2);
        h1i = 0.5*(*datai21-*dataj21);
        h2i = c2*(*datai2-*dataj2);
        h2r = -c2*(*datai21+*dataj21);
        *datai2 = h1r+wr*h2r-wi*h2i;
        *datai21 = h1i+wr*h2i+wi*h2r;
        *dataj2 = h1r-wr*h2r+wi*h2i;
        *dataj21 = -h1i+wr*h2i+wi*h2r;
        }
      }
    wr = (wtemp=wr)*wpr-wi*wpi+wr;
    wi = wi*wpr+wtemp*wpi+wi;
    }

  if (isign == -1)
    fourn(data, nn, 2, isign);

  return;
  }

