/*
*				fft.c
*
* Functions dealing with FFT.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 2008-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	PSFEx is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
* 	(at your option) any later version.
*	PSFEx is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with PSFEx.  If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		10/07/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef FFTW3_H
#include FFTW_H
#endif

#include "define.h"
#include "globals.h"
#include "fft.h"
#ifdef USE_THREADS
#include "threads.h"
#endif

 int    firsttimeflag;
#ifdef USE_THREADS
pthread_mutex_t	fftmutex;
#endif

/****** fft_init ************************************************************
PROTO	void fft_init(int nthreads)
PURPOSE	Initialize the FFT routines
INPUT	-.
OUTPUT	-.
NOTES	Global preferences are used for multhreading.
AUTHOR	E. Bertin (IAP)
VERSION	26/06/2009
 ***/
void    fft_init(int nthreads)
 {
  if (!firsttimeflag)
    {
#ifdef USE_THREADS
    QPTHREAD_MUTEX_INIT(&fftmutex, NULL);
#ifdef HAVE_FFTWF_MP
    if (nthreads > 1)
      {
      if (!fftwf_init_threads())
        error(EXIT_FAILURE, "*Error*: thread initialization failed in ","FFTW");
      fftwf_plan_with_nthreads(nthreads);
      }
#endif
#endif
    firsttimeflag = 1;
    }

  return;
  }


/****** fft_end ************************************************************
PROTO	void fft_init(int nthreads)
PURPOSE	Clear up stuff set by FFT routines
INPUT	-.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	26/06/2009
 ***/
void    fft_end(int nthreads)
 {

  if (firsttimeflag)
    {
    firsttimeflag = 0;
#ifdef USE_THREADS
    if (nthreads > 1)
      {
#ifdef HAVE_FFTWF_MP
      fftwf_cleanup_threads();
#endif
      QPTHREAD_MUTEX_DESTROY(&fftmutex);
      }
    else
#endif
      fftwf_cleanup();
    }

  return;
  }


/****** fft_conv ************************************************************
PROTO	void fft_conv(float *data1, float *fdata2, int width, int height)
PURPOSE	Optimized 2-dimensional FFT convolution using the FFTW library.
INPUT	ptr to the first image,
	ptr to the Fourier transform of the second image,
	image width,
	image height.
OUTPUT	-.
NOTES	For data1 and fdata2, memory must be allocated for
	size[0]* ... * 2*(size[naxis-1]/2+1) floats (padding required).
AUTHOR	E. Bertin (IAP)
VERSION	15/04/2009
 ***/
void    fft_conv(float *data1, float *fdata2, int width, int height)
  {
   fftwf_plan	plan;
   float	*fdata1,*fdata1p,*fdata2p,
		real,imag, fac;
   int		i, npix,npix2;

/* Convert axis indexing to that of FFTW */
  npix = width*height;
  npix2 = (((width>>1) + 1)<< 1) * height;

/* Forward FFT "in place" for data1 */
#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
  QFFTWMALLOC(fdata1, float, npix2);
  plan = fftwf_plan_dft_r2c_2d(height, width, data1,
        (fftwf_complex *)fdata1, FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif
  fftwf_execute(plan);

#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
  fftwf_destroy_plan(plan);
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif

/* Actual convolution (Fourier product) */
  fac = 1.0/npix;  
  fdata1p = fdata1;
  fdata2p = fdata2;
  for (i=npix2/2; i--; fdata2p+=2)
    {
    real = *fdata1p **fdata2p - *(fdata1p+1)**(fdata2p+1);
    imag = *(fdata1p+1)**fdata2p + *fdata1p**(fdata2p+1);
    *(fdata1p++) = fac*real;
    *(fdata1p++) = fac*imag;
    }

/* Reverse FFT */
#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
  plan = fftwf_plan_dft_c2r_2d(height, width, (fftwf_complex *)fdata1, 
        data1, FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif

  fftwf_execute(plan);

#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
  fftwf_destroy_plan(plan);
/* Free the fdata1 scratch array */
  QFFTWFREE(fdata1);
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif

  return;
  }


/****** fft_rtf ************************************************************
PROTO	float *fft_rtf(float *data, int width, int height)
PURPOSE	Optimized 2-dimensional FFT "in place" using the FFTW library.
INPUT	ptr to the image,
	image width,
	image height.
OUTPUT	Pointer to the compressed, memory-allocated Fourier transform.
NOTES	Input data may end up corrupted.
AUTHOR	E. Bertin (IAP)
VERSION	15/04/2009
 ***/
float	*fft_rtf(float *data, int width, int height)
  {
   fftwf_plan   plan;
   float	*fdata;
   int		npix2;

/* Convert axis indexing to that of FFTW */
  npix2 = (((width>>1) + 1)<< 1) * height;

/* Forward FFT "in place" for data1 */
#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
  QFFTWMALLOC(fdata, float, npix2);
  plan = fftwf_plan_dft_r2c_2d(height, width, data,
        (fftwf_complex *)fdata, FFTW_ESTIMATE);
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif

  fftwf_execute(plan);

#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
  fftwf_destroy_plan(plan);
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif

  return fdata;
  }



/****** fft_ctf ************************************************************
PROTO	void fft_ctf(float *data, int width, int height int sign)
PURPOSE	Optimized 2-dimensional complex FFT "in place" using the FFTW library.
INPUT	ptr to the image,
	image width,
	image height.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	15/04/2009
 ***/
void	fft_ctf(float *data, int width, int height, int sign)
  {
   fftwf_plan	plan;

/* Forward FFT "in place" for data1 */
#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
  plan = fftwf_plan_dft_2d(height, width, (fftwf_complex *)data,
		(fftwf_complex *)data, sign, FFTW_ESTIMATE);
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif

  fftwf_execute(plan);

#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
  fftwf_destroy_plan(plan);
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif

  return;
  }


/****** fft_shift ************************************************************
PROTO	void fft_shift(float *data, int width, int height)
PURPOSE	(de-)scramble 2-dimensional Fourier plane.
INPUT	ptr to the image,
	image width,
	image height.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	21/04/2009
 ***/
void	fft_shift(float *data, int width, int height)
  {
   float	*temp, *tempt, *datat;
   int		i, xc,yc, x,y,x2,y2, npix;

  npix = width*height;
  QMALLOC(temp, float, width*height);
  datat = data;
  xc = width/2;
  yc = height/2;
  for (y=0; y<height; y++)
    {
    y2 = y-yc;
    if (y2<0)
      y2 += height;
    y2 *= width;
    for (x=0; x<width; x++)
      {
      x2 = x-xc;
      if (x2<0)
        x2 += width;
      temp[x2+y2] = *(datat++);
      }
    }

  tempt = temp;
  for (i=npix; i--;)
    *(data++) = *(tempt++);

  free(temp);

  return;
  }



