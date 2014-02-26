/*
*				xml.c
*
* Handle XML metadata.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 2005-2014 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		26/02/2014
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "define.h"
#include "globals.h"
#include "fits/fitscat.h"
#include "cplot.h"
#include "key.h"
#include "field.h"
#include "prefs.h"
#include "psf.h"
#include "xml.h"

extern time_t		thetime,thetime2;	/* from makeit.c */
extern pkeystruct	key[];			/* from preflist.h */
extern char		keylist[][32];		/* from preflist.h */
 
fieldstruct		**field_xml;
int			nxml, nxmlmax;


/****** init_xml ************************************************************
PROTO	int init_xml(int ncat)
PURPOSE	Initialize a set of meta-data kept in memory before being written to the
	XML file
INPUT	Number of catalogues.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	04/07/2008
 ***/
int	init_xml(int ncat)
  {
  if (ncat)
    {
    QMALLOC(field_xml, fieldstruct *, ncat);
    }
  else
    field_xml = NULL;
  nxml = 0;
  nxmlmax = ncat;

  return EXIT_SUCCESS;
  }


/****** end_xml ************************************************************
PROTO	void end_xml(void)
PURPOSE	Free the set of meta-data kept in memory.
INPUT	-.
OUTPUT	.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	04/07/2008
 ***/
void	end_xml(void)
  {
  free(field_xml);

  return;
  }


/****** update_xml ***********************************************************
PROTO	int update_xml(fieldstruct *field)
PURPOSE	Update a set of meta-data kept in memory before being written to the
	XML file
INPUT	Pointer to the current field.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	Global preferences are used.
AUTHOR	E. Bertin (IAP)
VERSION	04/07/2008
 ***/
int	update_xml(fieldstruct *field)
  {
  field_xml[nxml] = field;
  nxml++;

  return EXIT_SUCCESS;
  }


/****** write_xml ************************************************************
PROTO	int	write_xml(char *filename)
PURPOSE	Save meta-data to an XML file/stream.
INPUT	XML file name.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	06/10/2006
 ***/
int	write_xml(char *filename)
  {
   FILE		*file;

  if (!(file = fopen(prefs.xml_name, "w")))
    return RETURN_ERROR;

  write_xml_header(file);
  write_xml_meta(file, (char *)NULL);

  fprintf(file, "</RESOURCE>\n");
  fprintf(file, "</VOTABLE>\n");

  fclose(file);

  return RETURN_OK;
  }


/****** write_xml_header ******************************************************
PROTO	int	write_xml_header(FILE *file)
PURPOSE	Save an XML-VOtable header to an XML file/stream
INPUT	file or stream pointer.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	Global preferences are used.
AUTHOR	E. Bertin (IAP)
VERSION	17/06/2012
 ***/
int	write_xml_header(FILE *file)
  {
   char		sysname[16];

  fprintf(file, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
  fprintf(file, "<?xml-stylesheet type=\"text/xsl\" href=\"%s\"?>\n",
	prefs.xsl_name);
  fprintf(file, "<VOTABLE version=\"1.1\"\n"
	" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
	" xsi:noNamespaceSchemaLocation=\"http://www.ivoa.net/xml/VOTable/v1.1\">\n");
  fprintf(file, "<DESCRIPTION>produced by %s</DESCRIPTION>\n", BANNER);
  fprintf(file, "<!-- VOTable description at "
	"http://www.ivoa.net/Documents/latest/VOT.html -->\n");
  fprintf(file, "<RESOURCE ID=\"%s\" name=\"%s\">\n", BANNER, BANNER);
  fprintf(file, " <DESCRIPTION>Data related to %s"
	"</DESCRIPTION>\n", BANNER);
  fprintf(file, " <INFO name=\"QUERY_STATUS\" value=\"OK\" />\n");
  sprintf(sysname, "ICRS");

  fprintf(file, " <COOSYS ID=\"J2000\" equinox=\"J2000\""
	" epoch=\"2000.0\" system=\"ICRS\"/>\n");

  return RETURN_OK;
  }


/****** write_xml_meta ********************************************************
PROTO	int	write_xml_meta(FILE *file, char *error)
PURPOSE	Save meta-data to an XML-VOTable file or stream
INPUT	Pointer to the output file (or stream),
	Pointer to an error msg (or NULL).
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	26/06/2012
 ***/
int	write_xml_meta(FILE *file, char *error)
  {
   fieldstruct		*field;
   psfstruct		*psf;
   struct tm		*tm;
   char			*pspath,*psuser, *pshost, *str;
   double		minrad_min,minrad_mean,minrad_max,
			sampling_min,sampling_mean,sampling_max,
			chi2_min,chi2_mean,chi2_max,
			fwhm_min,fwhm_mean,fwhm_max,
			fwhm_wcs_min,fwhm_wcs_mean,fwhm_wcs_max,
			ellipticity_min,ellipticity_mean,ellipticity_max,
			ellipticity1_min,ellipticity1_mean,ellipticity1_max,
			ellipticity2_min,ellipticity2_mean,ellipticity2_max,
			beta_min,beta_mean,beta_max,
			residuals_min,residuals_mean,residuals_max,
			pffwhm_min,pffwhm_mean,pffwhm_max,
			pffwhm_wcs_min,pffwhm_wcs_mean,pffwhm_wcs_max,
			pfellipticity_min,pfellipticity_mean,pfellipticity_max,
			pfellipticity1_min,pfellipticity1_mean,pfellipticity1_max,
			pfellipticity2_min,pfellipticity2_mean,pfellipticity2_max,
			pfbeta_min,pfbeta_mean,pfbeta_max,
			pfresiduals_min,pfresiduals_mean,pfresiduals_max,
			symresiduals_min,symresiduals_mean,symresiduals_max,
			noiseqarea_min,noiseqarea_mean,noiseqarea_max,
			pixscale_wcs_min,pixscale_wcs_mean,pixscale_wcs_max,
			nloaded_mean,naccepted_mean;
   int			d,n,e,
			nloaded_min,nloaded_max,nloaded_total,
			naccepted_min,naccepted_max,naccepted_total, neff, next;
#ifdef HAVE_PLPLOT
   char			plotfilename[MAXCHAR],
			*pstr;
   int			cp[CPLOT_NTYPES],
			j,t, nplot, pnplot, pngindex, pngflag;
#endif

/* Processing date and time if msg error present */
  if (error)
    {
    thetime2 = time(NULL);
    tm = localtime(&thetime2);
    sprintf(prefs.sdate_end,"%04d-%02d-%02d",
        tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
    sprintf(prefs.stime_end,"%02d:%02d:%02d",
        tm->tm_hour, tm->tm_min, tm->tm_sec);
    prefs.time_diff = difftime(thetime2, thetime);
    }

/* Username */
  psuser = pspath = pshost = NULL;
#ifdef HAVE_GETENV
  if (!(psuser=getenv("USERNAME")))	/* Cygwin,... */
    psuser = getenv("LOGNAME");		/* Linux,... */
  pspath = getenv("PWD");
  pshost = getenv("HOSTNAME");
#endif

  fprintf(file, " <RESOURCE ID=\"MetaData\" name=\"MetaData\">\n");
  fprintf(file, "  <DESCRIPTION>%s meta-data</DESCRIPTION>\n", BANNER);
  fprintf(file, "  <INFO name=\"QUERY_STATUS\" value=\"OK\" />\n");
  fprintf(file, "  <PARAM name=\"Software\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.title;meta.software\" value=\"%s\"/>\n",
	BANNER);
  fprintf(file, "  <PARAM name=\"Version\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.version;meta.software\" value=\"%s\"/>\n",
	MYVERSION);
  fprintf(file, "  <PARAM name=\"Soft_URL\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.ref.url;meta.software\" value=\"%s\"/>\n",
	WEBSITE);
  fprintf(file, "  <PARAM name=\"Soft_Auth\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.bib.author;meta.software\" value=\"%s\"/>\n",
	"Emmanuel Bertin");
  fprintf(file, "  <PARAM name=\"Soft_Ref\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.bib.bibcode;meta.software\" value=\"%s\"/>\n",
	"2006ASPC..351..112B");
  fprintf(file, "  <PARAM name=\"NThreads\" datatype=\"int\""
	" ucd=\"meta.number;meta.software\" value=\"%d\"/>\n",
    	prefs.nthreads);
  fprintf(file, "  <PARAM name=\"Date\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"time.event.end;meta.software\" value=\"%s\"/>\n",
	prefs.sdate_end);
  fprintf(file, "  <PARAM name=\"Time\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"time.event.end;meta.software\" value=\"%s\"/>\n",
	prefs.stime_end);
  fprintf(file, "  <PARAM name=\"Duration\" datatype=\"float\""
	" ucd=\"time.event;meta.software\" value=\"%.0f\" unit=\"s\"/>\n",
	prefs.time_diff);

  fprintf(file, "  <PARAM name=\"User\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.curation\" value=\"%s\"/>\n",
	psuser);
  fprintf(file, "  <PARAM name=\"Host\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.curation\" value=\"%s\"/>\n",
	pshost);
  fprintf(file, "  <PARAM name=\"Path\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.dataset\" value=\"%s\"/>\n",
	pspath);

  if (error)
    {
    fprintf(file, "\n  <!-- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	"!!!!!!!!!!!!!!!!!!!! -->\n");
    fprintf(file, "  <!-- !!!!!!!!!!!!!!!!!!!!!! an Error occured"
	" !!!!!!!!!!!!!!!!!!!!! -->\n");
    fprintf(file, "  <!-- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	"!!!!!!!!!!!!!!!!!!!! -->\n");
    fprintf(file,"  <PARAM name=\"Error_Msg\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta\" value=\"%s\"/>\n", error);
    fprintf(file, "  <!-- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	"!!!!!!!!!!!!!!!!!!!! -->\n");
    fprintf(file, "  <!-- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	"!!!!!!!!!!!!!!!!!!!! -->\n\n");
    }

/* Test if PNG plots are being produced */
#ifdef HAVE_PLPLOT
  nplot = pnplot = pngflag = 0;
  for (j=0; j<prefs.ncplot_device; j++)
    if ((prefs.cplot_device[j] == CPLOT_PNG))
      {
      pngflag = 1;
      break;
      }
#endif

/* PSF meta-data per field */
  fprintf(file, "  <TABLE ID=\"PSF_Fields\" name=\"PSF_Fields\">\n");
  fprintf(file, "   <DESCRIPTION>PSF metadata and stats per field gathered by "
	"%s</DESCRIPTION>\n", BANNER);
  fprintf(file, "   <!-- NFields may be 0"
	" if an error occurred early in the processing -->\n");
  fprintf(file, "   <PARAM name=\"NFields\" datatype=\"int\""
	" ucd=\"meta.number;meta.dataset\" value=\"%d\"/>\n", nxmlmax);
  fprintf(file, "   <!-- CurrField may differ from NFields"
	" if an error occurred -->\n");
  fprintf(file, "   <PARAM name=\"CurrField\" datatype=\"int\""
	" ucd=\"meta.number;meta.dataset\" value=\"%d\"/>\n", nxml);
  fprintf(file, "   <PARAM name=\"NSnapshots\" datatype=\"int\""
	" arraysize=\"%d\" ucd=\"meta.number;meta.dataset\""
	" value=\"%d",
	prefs.ncontext_group? prefs.ncontext_group : 1,
	prefs.ncontext_group? prefs.context_nsnap : 1);
  for (d=1; d<prefs.ncontext_group; d++)
    fprintf(file, " %d", prefs.context_nsnap);
  fprintf(file, "\"/>\n");

  fprintf(file, "   <FIELD name=\"Catalog_Name\" datatype=\"char\""
	" arraysize=\"*\" ucd=\"meta.id;meta.table;meta.file\"/>\n");
  fprintf(file, "   <FIELD name=\"Image_Ident\" datatype=\"char\""
	" arraysize=\"*\" ucd=\"meta.id;obs.field\"/>\n");
  fprintf(file, "   <FIELD name=\"NExtensions\" datatype=\"int\""
        " ucd=\"meta.number\"/>\n");
  fprintf(file, "   <FIELD name=\"NStars_Loaded_Total\" datatype=\"int\""
	" ucd=\"meta.number;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"NStars_Loaded_Min\" datatype=\"int\""
	" ucd=\"meta.number;stat.min;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"NStars_Loaded_Mean\" datatype=\"float\""
	" ucd=\"meta.number;stat.mean;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"NStars_Loaded_Max\" datatype=\"int\""
	" ucd=\"meta.number;stat.max;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"NStars_Accepted_Total\" datatype=\"int\""
	" ucd=\"meta.number;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"NStars_Accepted_Min\" datatype=\"int\""
	" ucd=\"meta.number;stat.min;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"NStars_Accepted_Mean\" datatype=\"int\""
	" ucd=\"meta.number;stat.mean;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"NStars_Accepted_Max\" datatype=\"int\""
	" ucd=\"meta.number;stat.max;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_FromFluxRadius_Min\" unit=\"pix\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_FromFluxRadius_Mean\" unit=\"pix\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_FromFluxRadius_Max\" unit=\"pix\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Sampling_Min\" unit=\"pix\""
	" datatype=\"float\""
	" ucd=\"arith.factor;instr.pixel;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Sampling_Mean\" unit=\"pix\""
	" datatype=\"float\""
	" ucd=\"arith.factor;instr.pixel;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Sampling_Max\" unit=\"pix\""
	" datatype=\"float\""
	" ucd=\"arith.factor;instr.pixel;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Chi2_Min\" datatype=\"float\""
	" ucd=\"stat.fit.chi2;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Chi2_Mean\" datatype=\"float\""
	" ucd=\"stat.fit.chi2;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Chi2_Max\" datatype=\"float\""
	" ucd=\"stat.fit.chi2;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_Min\" unit=\"pix\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_Mean\" unit=\"pix\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_Max\" unit=\"pix\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_WCS_Min\" unit=\"arcsec\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_WCS_Mean\" unit=\"arcsec\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_WCS_Max\" unit=\"arcsec\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity_Min\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity_Mean\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity_Max\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity1_Min\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity1_Mean\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity1_Max\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity2_Min\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity2_Mean\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity2_Max\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"MoffatBeta_Min\" datatype=\"float\""
	" ucd=\"stat.param;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"MoffatBeta_Mean\" datatype=\"float\""
	" ucd=\"stat.param;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"MoffatBeta_Max\" datatype=\"float\""
	" ucd=\"stat.param;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Residuals_Min\" datatype=\"float\""
	" ucd=\"stat.fit.residual;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Residuals_Mean\" datatype=\"float\""
	" ucd=\"stat.fit.residual;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Residuals_Max\" datatype=\"float\""
	" ucd=\"stat.fit.residual;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_PixelFree_Min\" unit=\"pix\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_PixelFree_Mean\" unit=\"pix\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_PixelFree_Max\" unit=\"pix\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_PixelFree_WCS_Min\" unit=\"arcsec\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_PixelFree_WCS_Mean\" unit=\"arcsec\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_PixelFree_WCS_Max\" unit=\"arcsec\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity_PixelFree_Min\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity_PixelFree_Mean\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity_PixelFree_Max\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity1_PixelFree_Min\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity1_PixelFree_Mean\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity1_PixelFree_Max\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity2_PixelFree_Min\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity2_PixelFree_Mean\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity2_PixelFree_Max\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"MoffatBeta_PixelFree_Min\" datatype=\"float\""
	" ucd=\"stat.param;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"MoffatBeta_PixelFree_Mean\" datatype=\"float\""
	" ucd=\"stat.param;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"MoffatBeta_PixelFree_Max\" datatype=\"float\""
	" ucd=\"stat.param;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Residuals_PixelFree_Min\" datatype=\"float\""
	" ucd=\"stat.fit.residual;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Residuals_PixelFree_Mean\" datatype=\"float\""
	" ucd=\"stat.fit.residual;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Residuals_PixelFree_Max\" datatype=\"float\""
	" ucd=\"stat.fit.residual;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Asymmetry_Min\" datatype=\"float\""
	" ucd=\"stat.fit.residual;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Asymmetry_Mean\" datatype=\"float\""
	" ucd=\"stat.fit.residual;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Asymmetry_Max\" datatype=\"float\""
	" ucd=\"stat.fit.residual;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Area_Noise_Min\" datatype=\"float\""
	" ucd=\"phys.area;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Area_Noise_Mean\" datatype=\"float\""
	" ucd=\"phys.area;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Area_Noise_Max\" datatype=\"float\""
	" ucd=\"phys.area;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"PixelScale_WCS_Min\" unit=\"arcsec\""
	" datatype=\"float\""
	" ucd=\"instr.scale;stat.min;instr.pixel\"/>\n");
  fprintf(file, "   <FIELD name=\"PixelScale_WCS_Mean\" unit=\"arcsec\""
	" datatype=\"float\""
	" ucd=\"instr.scale;stat.mean;instr.pixel\"/>\n");
  fprintf(file, "   <FIELD name=\"PixelScale_WCS_Max\" unit=\"arcsec\""
	" datatype=\"float\""
	" ucd=\"instr.scale;stat.max;instr.pixel\"/>\n");

/*-- Checkplots */
#ifdef HAVE_PLPLOT
  if (pngflag)
    {
    pnplot = nplot;
    if ((pngindex=cplot_check(CPLOT_COUNTS)) != RETURN_ERROR)
      {
      fprintf(file, "   <FIELD name=\"Plot_Counts\" datatype=\"char\""
        " arraysize=\"*\" ucd=\"meta.id;meta.file\"/>\n");
      cp[nplot++] = pngindex;
      }
    if ((pngindex=cplot_check(CPLOT_COUNTFRAC)) != RETURN_ERROR)
      {
      fprintf(file, "   <FIELD name=\"Plot_Count_Fraction\" datatype=\"char\""
        " arraysize=\"*\" ucd=\"meta.id;meta.file\"/>\n");
      cp[nplot++] = pngindex;
      }
    if ((pngindex=cplot_check(CPLOT_FWHM)) != RETURN_ERROR)
      {
      fprintf(file, "   <FIELD name=\"Plot_FWHM\" datatype=\"char\""
        " arraysize=\"*\" ucd=\"meta.id;meta.file\"/>\n");
      cp[nplot++] = pngindex;
      }
    if ((pngindex=cplot_check(CPLOT_ELLIPTICITY)) != RETURN_ERROR)
      {
      fprintf(file, "   <FIELD name=\"Plot_Ellipticity\" datatype=\"char\""
        " arraysize=\"*\" ucd=\"meta.id;meta.file\"/>\n");
      cp[nplot++] = pngindex;
      }
    }
#endif

  fprintf(file, "   <DATA><TABLEDATA>\n");
  for (n=0; n<nxml; n++)
    {
/*-- Compute min,average and max of Moffat fitted parameters */
    nloaded_min = naccepted_min = 2<<29;
    nloaded_max = naccepted_max = nloaded_total = naccepted_total = 0;
    minrad_min = sampling_min = chi2_min = fwhm_min = fwhm_wcs_min
	= ellipticity_min = ellipticity1_min = ellipticity2_min
	= beta_min = residuals_min = pffwhm_min = pffwhm_wcs_min
	= pfellipticity_min = pfellipticity1_min = pfellipticity2_min
	= pfbeta_min = pfresiduals_min = symresiduals_min
	= noiseqarea_min = pixscale_wcs_min = BIG;
    minrad_mean = sampling_mean = chi2_mean = fwhm_mean = fwhm_wcs_mean
	= ellipticity_mean = ellipticity1_mean = ellipticity2_mean
	= beta_mean = residuals_mean = pffwhm_mean = pffwhm_wcs_mean
	= pfellipticity_mean = pfellipticity1_mean = pfellipticity2_mean
	= pfbeta_mean = pfresiduals_mean = symresiduals_mean = noiseqarea_mean
	= pixscale_wcs_mean = nloaded_mean = naccepted_mean = 0.0;
    minrad_max = sampling_max = chi2_max = fwhm_max = fwhm_wcs_max
	= ellipticity_max = ellipticity1_max = ellipticity2_max
	= beta_max = residuals_max = pffwhm_max = pffwhm_wcs_max
	= pfellipticity_max = pfellipticity1_max = pfellipticity2_max
	= pfbeta_max = pfresiduals_max = symresiduals_max
	= noiseqarea_max = pixscale_wcs_max = -BIG;
    neff = 0;
    field = field_xml[n];
    for (e=0; e<field->next; e++)
      {
      psf = field->psf[e];
      nloaded_total += psf->samples_loaded;
      if (psf->samples_loaded < nloaded_min)
        nloaded_min = psf->samples_loaded ;
      nloaded_mean += (double)psf->samples_loaded;
      if (psf->samples_loaded > nloaded_max)
        nloaded_max = psf->samples_loaded;
      naccepted_total += psf->samples_accepted;
      if (psf->samples_accepted < naccepted_min)
        naccepted_min = psf->samples_accepted;
      naccepted_mean += (double)psf->samples_accepted;
      if (psf->samples_accepted > naccepted_max)
        naccepted_max = psf->samples_accepted ;
/*---- Drop it if no valid stars have been kept */
      if (!psf->samples_accepted)
        continue;
      neff++;
      if (psf->fwhm < minrad_min)
        minrad_min = psf->fwhm;
      minrad_mean += psf->fwhm;
      if (psf->fwhm > minrad_max)
        minrad_max = psf->fwhm;
      if (psf->pixstep < sampling_min)
        sampling_min = psf->pixstep;
      sampling_mean += psf->pixstep;
      if (psf->pixstep > sampling_max)
        sampling_max = psf->pixstep;
      if (psf->chi2 < chi2_min)
        chi2_min = psf->chi2;
      chi2_mean += psf->chi2;
      if (psf->chi2 > chi2_max)
        chi2_max = psf->chi2;
/*---- Moffat fit */
      if (psf->moffat_fwhm_min < fwhm_min)
        fwhm_min = psf->moffat_fwhm_min;
      fwhm_mean += psf->moffat_fwhm;
      if (psf->moffat_fwhm_max > fwhm_max)
        fwhm_max = psf->moffat_fwhm_max;
      if (psf->moffat_fwhm_wcs_min < fwhm_wcs_min)
        fwhm_wcs_min = psf->moffat_fwhm_wcs_min;
      fwhm_wcs_mean += psf->moffat_fwhm_wcs;
      if (psf->moffat_fwhm_wcs_max > fwhm_wcs_max)
        fwhm_wcs_max = psf->moffat_fwhm_wcs_max;
      if (psf->moffat_ellipticity_min < ellipticity_min)
        ellipticity_min = psf->moffat_ellipticity_min;
      ellipticity_mean += psf->moffat_ellipticity;
      if (psf->moffat_ellipticity_max > ellipticity_max)
        ellipticity_max = psf->moffat_ellipticity_max;
      if (psf->moffat_ellipticity1_min < ellipticity1_min)
        ellipticity1_min = psf->moffat_ellipticity1_min;
      ellipticity1_mean += psf->moffat_ellipticity1;
      if (psf->moffat_ellipticity1_max > ellipticity1_max)
        ellipticity1_max = psf->moffat_ellipticity1_max;
      if (psf->moffat_ellipticity2_min < ellipticity2_min)
        ellipticity2_min = psf->moffat_ellipticity2_min;
      ellipticity2_mean += psf->moffat_ellipticity2;
      if (psf->moffat_ellipticity2_max > ellipticity2_max)
        ellipticity2_max = psf->moffat_ellipticity2_max;
      if (psf->moffat_beta_min < beta_min)
        beta_min = psf->moffat_beta_min;
      beta_mean += psf->moffat_beta;
      if (psf->moffat_beta_max > beta_max)
        beta_max = psf->moffat_beta_max;
      if (psf->moffat_residuals_min < residuals_min)
        residuals_min = psf->moffat_residuals_min;
      residuals_mean += psf->moffat_residuals;
      if (psf->moffat_residuals_max > residuals_max)
        residuals_max = psf->moffat_residuals_max;
/*---- Pixel-free Moffat fit */
      if (psf->pfmoffat_fwhm_min < pffwhm_min)
        pffwhm_min = psf->pfmoffat_fwhm_min;
      pffwhm_mean += psf->pfmoffat_fwhm;
      if (psf->pfmoffat_fwhm_max > pffwhm_max)
        pffwhm_max = psf->pfmoffat_fwhm_max;
      if (psf->pfmoffat_fwhm_wcs_min < pffwhm_wcs_min)
        pffwhm_wcs_min = psf->pfmoffat_fwhm_wcs_min;
      pffwhm_wcs_mean += psf->pfmoffat_fwhm_wcs;
      if (psf->pfmoffat_fwhm_wcs_max > pffwhm_wcs_max)
        pffwhm_wcs_max = psf->pfmoffat_fwhm_wcs_max;
      if (psf->pfmoffat_ellipticity_min < pfellipticity_min)
        pfellipticity_min = psf->pfmoffat_ellipticity_min;
      pfellipticity_mean += psf->pfmoffat_ellipticity;
      if (psf->pfmoffat_ellipticity_max > pfellipticity_max)
        pfellipticity_max = psf->pfmoffat_ellipticity_max;
      if (psf->pfmoffat_ellipticity1_min < pfellipticity1_min)
        pfellipticity1_min = psf->pfmoffat_ellipticity1_min;
      pfellipticity1_mean += psf->pfmoffat_ellipticity1;
      if (psf->pfmoffat_ellipticity1_max > pfellipticity1_max)
        pfellipticity1_max = psf->pfmoffat_ellipticity1_max;
      if (psf->pfmoffat_ellipticity2_min < pfellipticity2_min)
        pfellipticity2_min = psf->pfmoffat_ellipticity2_min;
      pfellipticity2_mean += psf->pfmoffat_ellipticity2;
      if (psf->pfmoffat_ellipticity2_max > pfellipticity2_max)
        pfellipticity2_max = psf->pfmoffat_ellipticity2_max;
      if (psf->pfmoffat_beta_min < pfbeta_min)
        pfbeta_min = psf->pfmoffat_beta_min;
      pfbeta_mean += psf->pfmoffat_beta;
      if (psf->pfmoffat_beta_max > pfbeta_max)
        pfbeta_max = psf->pfmoffat_beta_max;
      if (psf->pfmoffat_residuals_min < pfresiduals_min)
        pfresiduals_min = psf->pfmoffat_residuals_min;
      pfresiduals_mean += psf->pfmoffat_residuals;
      if (psf->pfmoffat_residuals_max > pfresiduals_max)
        pfresiduals_max = psf->pfmoffat_residuals_max;
/*---- Asymmetry measurement */
      if (psf->sym_residuals_min < symresiduals_min)
        symresiduals_min = psf->sym_residuals_min;
      symresiduals_mean += psf->sym_residuals;
      if (psf->sym_residuals_max > symresiduals_max)
        symresiduals_max = psf->sym_residuals_max;
/*---- Noise equal area measurement */
      if (psf->noiseqarea_min < noiseqarea_min)
        noiseqarea_min = psf->noiseqarea_min;
      noiseqarea_mean += psf->noiseqarea;
      if (psf->noiseqarea_max > noiseqarea_max)
        noiseqarea_max = psf->noiseqarea_max;
      if (psf->pixscale_wcs_min < pixscale_wcs_min)
        pixscale_wcs_min = psf->pixscale_wcs_min;
      pixscale_wcs_mean += psf->pixscale_wcs;
      if (psf->pixscale_wcs_max > pixscale_wcs_max)
        pixscale_wcs_max = psf->pixscale_wcs_max;
      }

    if (field->next>1)
      {
      nloaded_mean /= (double)field->next;
      naccepted_mean /= (double)field->next;
      }

    if (neff>1)
      {
      minrad_mean /= (double)neff;
      sampling_mean /= (double)neff;
      chi2_mean /= (double)neff;
      fwhm_mean /= (double)neff;
      fwhm_wcs_mean /= (double)neff;
      ellipticity_mean /= (double)neff;
      ellipticity1_mean /= (double)neff;
      ellipticity2_mean /= (double)neff;
      beta_mean /= (double)neff;
      residuals_mean /= (double)neff;
      pffwhm_mean /= (double)neff;
      pffwhm_wcs_mean /= (double)neff;
      pfellipticity_mean /= (double)neff;
      pfellipticity1_mean /= (double)neff;
      pfellipticity2_mean /= (double)neff;
      pfbeta_mean /= (double)neff;
      pfresiduals_mean /= (double)neff;
      symresiduals_mean /= (double)neff;
      noiseqarea_mean /= (double)neff;
      pixscale_wcs_mean /= (double)neff;
      }
    else if (neff==0)
      minrad_min = sampling_min = chi2_min = fwhm_min = fwhm_wcs_min
	= ellipticity_min = ellipticity1_min = ellipticity2_min
	= beta_min = residuals_min = pffwhm_min = pffwhm_wcs_min
	= pfellipticity_min = pfellipticity1_min = pfellipticity2_min
	= pfbeta_min = pfresiduals_min = symresiduals_min = noiseqarea_min
	= pixscale_wcs_min
	= minrad_mean = sampling_mean = chi2_mean = fwhm_mean = fwhm_wcs_mean
	= ellipticity_mean = ellipticity1_mean = ellipticity2_mean
	= beta_mean = residuals_mean = pffwhm_mean = pffwhm_wcs_mean
	= pfellipticity_mean = pfellipticity1_mean = pfellipticity2_mean
	= pfbeta_mean = pfresiduals_mean = symresiduals_mean = noiseqarea_mean
	= pixscale_wcs_mean
	= minrad_max = sampling_max = chi2_max = fwhm_max = fwhm_wcs_max
	= ellipticity_max  = ellipticity1_max = ellipticity2_max
	= beta_max = residuals_max = pffwhm_max = pffwhm_wcs_max
	= pfellipticity_max = pfellipticity1_max = pfellipticity2_max
	= pfbeta_max = pfresiduals_max = symresiduals_max = noiseqarea_max
	= pixscale_wcs_max = 0.0;

    fprintf(file, "    <TR>\n"
	"     <TD>%s</TD><TD>%s</TD><TD>%d</TD>\n"
        "     <TD>%d</TD><TD>%d</TD><TD>%.6g</TD><TD>%d</TD>\n"
        "     <TD>%d</TD><TD>%d</TD><TD>%.6g</TD><TD>%d</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n",
	field->rcatname, field->ident, field->next,
	nloaded_total, nloaded_min, nloaded_mean, nloaded_max,
	naccepted_total, naccepted_min, naccepted_mean, naccepted_max,
	minrad_min, minrad_mean, minrad_max,
	sampling_min, sampling_mean, sampling_max,
	chi2_min, chi2_mean, chi2_max,
	fwhm_min, fwhm_mean, fwhm_max,
	fwhm_wcs_min, fwhm_wcs_mean, fwhm_wcs_max,
	ellipticity_min, ellipticity_mean, ellipticity_max,
	ellipticity1_min, ellipticity1_mean, ellipticity1_max,
	ellipticity2_min, ellipticity2_mean, ellipticity2_max,
	beta_min, beta_mean, beta_max,
	residuals_min, residuals_mean, residuals_max,
	pffwhm_min, pffwhm_mean, pffwhm_max,
	pffwhm_wcs_min, pffwhm_wcs_mean, pffwhm_wcs_max,
	pfellipticity_min, pfellipticity_mean, pfellipticity_max,
	pfellipticity1_min, pfellipticity1_mean, pfellipticity1_max,
	pfellipticity2_min, pfellipticity2_mean, pfellipticity2_max,
	pfbeta_min, pfbeta_mean, pfbeta_max,
	pfresiduals_min, pfresiduals_mean, pfresiduals_max,
	symresiduals_min, symresiduals_mean, symresiduals_max,
	noiseqarea_min, noiseqarea_mean, noiseqarea_max,
	pixscale_wcs_min, pixscale_wcs_mean, pixscale_wcs_max);

/*-- Check-plots */
#ifdef HAVE_PLPLOT
    if (pngflag)
      {
      for (t=pnplot; t<nplot; t++)
        {
        strcpy(plotfilename, field->rcatname);
        if (!(pstr = strrchr(plotfilename, '.')))
          pstr = plotfilename+strlen(plotfilename);
        sprintf(pstr, ".png");
        fprintf(file, "     <TD>%s_%s</TD>\n",
		prefs.cplot_name[cp[t]], plotfilename);        
        }
      }
#endif
    fprintf(file, "    </TR>\n");
    }
  fprintf(file, "   </TABLEDATA></DATA>\n");
  fprintf(file, "  </TABLE>\n");

/* PSF meta-data per extension*/
  next = nxml? field_xml[0]->next : 0;
  fprintf(file, "  <TABLE ID=\"PSF_Extensions\" name=\"PSF_Extensions\">\n");
  fprintf(file, "   <DESCRIPTION>PSF metadata and stats per extension gathered"
	" by %s</DESCRIPTION>\n", BANNER);
  fprintf(file, "   <!-- NExtensions may be 0"
	" if an error occurred early in the processing -->\n");
  fprintf(file, "   <PARAM name=\"NExtensions\" datatype=\"int\""
	" ucd=\"meta.number;meta.dataset\" value=\"%d\"/>\n", next);
  fprintf(file, "   <FIELD name=\"Extension\" datatype=\"int\""
        " ucd=\"meta.record\"/>\n");
  fprintf(file, "   <FIELD name=\"NStars_Loaded_Total\" datatype=\"int\""
	" ucd=\"meta.number;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"NStars_Loaded_Min\" datatype=\"int\""
	" ucd=\"meta.number;stat.min;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"NStars_Loaded_Mean\" datatype=\"float\""
	" ucd=\"meta.number;stat.mean;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"NStars_Loaded_Max\" datatype=\"int\""
	" ucd=\"meta.number;stat.max;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"NStars_Accepted_Total\" datatype=\"int\""
	" ucd=\"meta.number;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"NStars_Accepted_Min\" datatype=\"int\""
	" ucd=\"meta.number;stat.min;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"NStars_Accepted_Mean\" datatype=\"float\""
	" ucd=\"meta.number;stat.mean;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"NStars_Accepted_Max\" datatype=\"int\""
	" ucd=\"meta.number;stat.max;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_FromFluxRadius_Min\" unit=\"pix\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_FromFluxRadius_Mean\" unit=\"pix\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_FromFluxRadius_Max\" unit=\"pix\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Sampling_Min\" unit=\"pix\""
	" datatype=\"float\""
	" ucd=\"arith.factor;instr.pixel;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Sampling_Mean\" unit=\"pix\""
	" datatype=\"float\""
	" ucd=\"arith.factor;instr.pixel;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Sampling_Max\" unit=\"pix\""
	" datatype=\"float\""
	" ucd=\"arith.factor;instr.pixel;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Chi2_Min\" datatype=\"float\""
	" ucd=\"stat.fit.chi2;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Chi2_Mean\" datatype=\"float\""
	" ucd=\"stat.fit.chi2;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Chi2_Max\" datatype=\"float\""
	" ucd=\"stat.fit.chi2;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_Min\" unit=\"pix\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_Mean\" unit=\"pix\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_Max\" unit=\"pix\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_WCS_Min\" unit=\"arcsec\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_WCS_Mean\" unit=\"arcsec\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_WCS_Max\" unit=\"arcsec\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity_Min\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity_Mean\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity_Max\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity1_Min\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity1_Mean\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity1_Max\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity2_Min\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity2_Mean\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity2_Max\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"MoffatBeta_Min\" datatype=\"float\""
	" ucd=\"stat.param;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"MoffatBeta_Mean\" datatype=\"float\""
	" ucd=\"stat.param;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"MoffatBeta_Max\" datatype=\"float\""
	" ucd=\"stat.param;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Residuals_Min\" datatype=\"float\""
	" ucd=\"stat.fit.residual;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Residuals_Mean\" datatype=\"float\""
	" ucd=\"stat.fit.residual;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Residuals_Max\" datatype=\"float\""
	" ucd=\"stat.fit.residual;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_PixelFree_Min\" unit=\"pix\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_PixelFree_Mean\" unit=\"pix\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_PixelFree_Max\" unit=\"pix\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_PixelFree_WCS_Min\" unit=\"arcsec\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_PixelFree_WCS_Mean\" unit=\"arcsec\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_PixelFree_WCS_Max\" unit=\"arcsec\""
	" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity_PixelFree_Min\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity_PixelFree_Mean\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity_PixelFree_Max\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity1_PixelFree_Min\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity1_PixelFree_Mean\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity1_PixelFree_Max\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity2_PixelFree_Min\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity2_PixelFree_Mean\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Ellipticity2_PixelFree_Max\" datatype=\"float\""
	" ucd=\"src.ellipticity;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"MoffatBeta_PixelFree_Min\" datatype=\"float\""
	" ucd=\"stat.param;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"MoffatBeta_PixelFree_Mean\" datatype=\"float\""
	" ucd=\"stat.param;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"MoffatBeta_PixelFree_Max\" datatype=\"float\""
	" ucd=\"stat.param;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Residuals_PixelFree_Min\" datatype=\"float\""
	" ucd=\"stat.fit.residual;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Residuals_PixelFree_Mean\" datatype=\"float\""
	" ucd=\"stat.fit.residual;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Residuals_PixelFree_Max\" datatype=\"float\""
	" ucd=\"stat.fit.residual;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Asymmetry_Min\" datatype=\"float\""
	" ucd=\"stat.fit.residual;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Asymmetry_Mean\" datatype=\"float\""
	" ucd=\"stat.fit.residual;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Asymmetry_Max\" datatype=\"float\""
	" ucd=\"stat.fit.residual;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Area_Noise_Min\" datatype=\"float\""
	" ucd=\"phys.area;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Area_Noise_Mean\" datatype=\"float\""
	" ucd=\"phys.area;stat.mean;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Area_Noise_Max\" datatype=\"float\""
	" ucd=\"phys.area;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"PixelScale_WCS_Min\" unit=\"arcsec\""
	" datatype=\"float\""
	" ucd=\"instr.scale;stat.min;instr.pixel\"/>\n");
  fprintf(file, "   <FIELD name=\"PixelScale_WCS_Mean\" unit=\"arcsec\""
	" datatype=\"float\""
	" ucd=\"instr.scale;stat.mean;instr.pixel\"/>\n");
  fprintf(file, "   <FIELD name=\"PixelScale_WCS_Max\" unit=\"arcsec\""
	" datatype=\"float\""
	" ucd=\"instr.scale;stat.max;instr.pixel\"/>\n");

  fprintf(file, "   <DATA><TABLEDATA>\n");
  for (e=0; e<next; e++)
    {
/*-- Compute min,average and max of Moffat fitted parameters */
    nloaded_min = naccepted_min = 2<<29;
    nloaded_max = naccepted_max = nloaded_total = naccepted_total = 0;
    minrad_min = sampling_min = chi2_min = fwhm_min = fwhm_wcs_min
	= ellipticity_min = ellipticity1_min = ellipticity2_min
	= beta_min = residuals_min = pffwhm_min
	= pfellipticity_min = pfellipticity1_min = pfellipticity2_min
	= pfbeta_min = pfresiduals_min = symresiduals_min
	= noiseqarea_min = pixscale_wcs_min = BIG;
    minrad_mean = sampling_mean = chi2_mean = fwhm_mean = fwhm_wcs_mean
	= ellipticity_mean = ellipticity1_mean = ellipticity2_mean
	= beta_mean = residuals_mean = pffwhm_mean = pffwhm_wcs_mean
	= pfellipticity_mean = pfellipticity1_mean = pfellipticity2_mean
	= pfbeta_mean = pfresiduals_mean = symresiduals_mean = noiseqarea_mean
	= pixscale_wcs_mean = nloaded_mean = naccepted_mean = 0.0;
    minrad_max = sampling_max = chi2_max = fwhm_max = fwhm_wcs_max
	= ellipticity_max = ellipticity1_max = ellipticity2_max
	= beta_max = residuals_max = pffwhm_max
	= pfellipticity_max = pfellipticity1_max = pfellipticity2_max
	= pfbeta_max = pfresiduals_max = symresiduals_max
	= noiseqarea_max = pixscale_wcs_max = -BIG;
    neff = 0;
    for (n=0; n<nxml; n++)
      {
      field = field_xml[n];
      psf = field->psf[e];
      nloaded_total += psf->samples_loaded;
      if (psf->samples_loaded < nloaded_min)
        nloaded_min = psf->samples_loaded ;
      nloaded_mean += (double)psf->samples_loaded;
      if (psf->samples_loaded > nloaded_max)
        nloaded_max = psf->samples_loaded;
      naccepted_total += psf->samples_accepted;
      if (psf->samples_accepted < naccepted_min)
        naccepted_min = psf->samples_accepted;
      naccepted_mean += (double)psf->samples_accepted;
      if (psf->samples_accepted > naccepted_max)
        naccepted_max = psf->samples_accepted ;
/*---- Drop it if no valid stars have been kept */
      if (!psf->samples_accepted)
        continue;
      neff++;
      if (psf->fwhm < minrad_min)
        minrad_min = psf->fwhm;
      minrad_mean += psf->fwhm;
      if (psf->fwhm > minrad_max)
        minrad_max = psf->fwhm;
      if (psf->pixstep < sampling_min)
        sampling_min = psf->pixstep;
      sampling_mean += psf->pixstep;
      if (psf->pixstep > sampling_max)
        sampling_max = psf->pixstep;
      if (psf->chi2 < chi2_min)
        chi2_min = psf->chi2;
      chi2_mean += psf->chi2;
      if (psf->chi2 > chi2_max)
        chi2_max = psf->chi2;
/*---- Moffat fit */
      if (psf->moffat_fwhm_min < fwhm_min)
        fwhm_min = psf->moffat_fwhm_min;
      fwhm_mean += psf->moffat_fwhm;
      if (psf->moffat_fwhm_max > fwhm_max)
        fwhm_max = psf->moffat_fwhm_max;
      if (psf->moffat_fwhm_wcs_min < fwhm_wcs_min)
        fwhm_wcs_min = psf->moffat_fwhm_wcs_min;
      fwhm_wcs_mean += psf->moffat_fwhm_wcs;
      if (psf->moffat_fwhm_wcs_max > fwhm_wcs_max)
        fwhm_wcs_max = psf->moffat_fwhm_wcs_max;
      if (psf->moffat_ellipticity_min < ellipticity_min)
        ellipticity_min = psf->moffat_ellipticity_min;
      ellipticity_mean += psf->moffat_ellipticity;
      if (psf->moffat_ellipticity_max > ellipticity_max)
        ellipticity_max = psf->moffat_ellipticity_max;
      if (psf->moffat_ellipticity1_min < ellipticity1_min)
        ellipticity1_min = psf->moffat_ellipticity1_min;
      ellipticity1_mean += psf->moffat_ellipticity1;
      if (psf->moffat_ellipticity1_max > ellipticity1_max)
        ellipticity1_max = psf->moffat_ellipticity1_max;
      if (psf->moffat_ellipticity2_min < ellipticity2_min)
        ellipticity2_min = psf->moffat_ellipticity2_min;
      ellipticity2_mean += psf->moffat_ellipticity2;
      if (psf->moffat_ellipticity2_max > ellipticity2_max)
        ellipticity2_max = psf->moffat_ellipticity2_max;
      if (psf->moffat_beta_min < beta_min)
        beta_min = psf->moffat_beta_min;
      beta_mean += psf->moffat_beta;
      if (psf->moffat_beta_max > beta_max)
        beta_max = psf->moffat_beta_max;
      if (psf->moffat_residuals_min < residuals_min)
        residuals_min = psf->moffat_residuals_min;
      residuals_mean += psf->moffat_residuals;
      if (psf->moffat_residuals_max > residuals_max)
        residuals_max = psf->moffat_residuals_max;
/*---- Pixel-free Moffat fit */
      if (psf->pfmoffat_fwhm_min < pffwhm_min)
        pffwhm_min = psf->pfmoffat_fwhm_min;
      pffwhm_mean += psf->pfmoffat_fwhm;
      if (psf->pfmoffat_fwhm_max > pffwhm_max)
        pffwhm_max = psf->pfmoffat_fwhm_max;
      if (psf->pfmoffat_fwhm_wcs_min < pffwhm_wcs_min)
        pffwhm_wcs_min = psf->pfmoffat_fwhm_wcs_min;
      pffwhm_wcs_mean += psf->pfmoffat_fwhm_wcs;
      if (psf->pfmoffat_fwhm_wcs_max > pffwhm_wcs_max)
        pffwhm_wcs_max = psf->pfmoffat_fwhm_wcs_max;
      if (psf->pfmoffat_ellipticity_min < pfellipticity_min)
        pfellipticity_min = psf->pfmoffat_ellipticity_min;
      pfellipticity_mean += psf->pfmoffat_ellipticity;
      if (psf->pfmoffat_ellipticity_max > pfellipticity_max)
        pfellipticity_max = psf->pfmoffat_ellipticity_max;
      if (psf->pfmoffat_ellipticity1_min < pfellipticity1_min)
        pfellipticity1_min = psf->pfmoffat_ellipticity1_min;
      pfellipticity1_mean += psf->pfmoffat_ellipticity1;
      if (psf->pfmoffat_ellipticity1_max > pfellipticity1_max)
        pfellipticity1_max = psf->pfmoffat_ellipticity1_max;
      if (psf->pfmoffat_ellipticity2_min < pfellipticity2_min)
        pfellipticity2_min = psf->pfmoffat_ellipticity2_min;
      pfellipticity2_mean += psf->pfmoffat_ellipticity2;
      if (psf->pfmoffat_ellipticity2_max > pfellipticity2_max)
        pfellipticity2_max = psf->pfmoffat_ellipticity2_max;
      if (psf->pfmoffat_beta_min < pfbeta_min)
        pfbeta_min = psf->pfmoffat_beta_min;
      pfbeta_mean += psf->pfmoffat_beta;
      if (psf->pfmoffat_beta_max > pfbeta_max)
        pfbeta_max = psf->pfmoffat_beta_max;
      if (psf->pfmoffat_residuals_min < pfresiduals_min)
        pfresiduals_min = psf->pfmoffat_residuals_min;
      pfresiduals_mean += psf->pfmoffat_residuals;
      if (psf->pfmoffat_residuals_max > pfresiduals_max)
        pfresiduals_max = psf->pfmoffat_residuals_max;
/*---- Asymmetry measurement */
      if (psf->sym_residuals_min < symresiduals_min)
        symresiduals_min = psf->sym_residuals_min;
      symresiduals_mean += psf->sym_residuals;
      if (psf->sym_residuals_max > symresiduals_max)
        symresiduals_max = psf->sym_residuals_max;
/*---- Noise equal area measurement */
      if (psf->noiseqarea_min < noiseqarea_min)
        noiseqarea_min = psf->noiseqarea_min;
      noiseqarea_mean += psf->noiseqarea;
      if (psf->noiseqarea_max > noiseqarea_max)
        noiseqarea_max = psf->noiseqarea_max;
      if (psf->pixscale_wcs_min < pixscale_wcs_min)
        pixscale_wcs_min = psf->pixscale_wcs_min;
      pixscale_wcs_mean += psf->pixscale_wcs;
      if (psf->pixscale_wcs_max > pixscale_wcs_max)
        pixscale_wcs_max = psf->pixscale_wcs_max;
      }

    if (nxml>1)
      {
      nloaded_mean /= (double)nxml;
      naccepted_mean /= (double)nxml;
      }

    if (neff>1)
      {
      minrad_mean /= (double)neff;
      sampling_mean /= (double)neff;
      chi2_mean /= (double)neff;
      fwhm_mean /= (double)neff;
      fwhm_wcs_mean /= (double)neff;
      ellipticity_mean /= (double)neff;
      ellipticity1_mean /= (double)neff;
      ellipticity2_mean /= (double)neff;
      beta_mean /= (double)neff;
      residuals_mean /= (double)neff;
      pffwhm_mean /= (double)neff;
      pffwhm_wcs_mean /= (double)neff;
      pfellipticity_mean /= (double)neff;
      pfellipticity1_mean /= (double)neff;
      pfellipticity2_mean /= (double)neff;
      pfbeta_mean /= (double)neff;
      pfresiduals_mean /= (double)neff;
      symresiduals_mean /= (double)neff;
      noiseqarea_mean /= (double)neff;
      pixscale_wcs_mean /= (double)neff;
      }
    else if (neff==0)
      minrad_min = sampling_min = chi2_min = fwhm_min = fwhm_wcs_min
	= ellipticity_min = ellipticity1_min = ellipticity2_min
	= beta_min = residuals_min = pffwhm_min = pffwhm_wcs_min
	= pfellipticity_min = pfellipticity1_min = pfellipticity2_min
	= pfbeta_min = pfresiduals_min = symresiduals_min = noiseqarea_min
	= pixscale_wcs_min
	= minrad_mean = sampling_mean = chi2_mean = fwhm_mean = fwhm_wcs_mean
	= ellipticity_mean = ellipticity1_mean = ellipticity2_mean
	= beta_mean = residuals_mean = pffwhm_mean = pffwhm_wcs_mean
	= pfellipticity_mean = pfellipticity1_mean = pfellipticity2_mean
	= pfbeta_mean = pfresiduals_mean = symresiduals_mean = noiseqarea_mean
	= pixscale_wcs_mean
	= minrad_max = sampling_max = chi2_max = fwhm_max = fwhm_wcs_max
	= ellipticity_max  = ellipticity1_max = ellipticity2_max
	= beta_max = residuals_max = pffwhm_max = pffwhm_wcs_max
	= pfellipticity_max = pfellipticity1_max = pfellipticity2_max
	= pfbeta_max = pfresiduals_max = symresiduals_max = noiseqarea_max
	= pixscale_wcs_max = 0.0;

    fprintf(file, "    <TR>\n"
	"     <TD>%d</TD>\n"
        "     <TD>%d</TD><TD>%d</TD><TD>%.6g</TD><TD>%d</TD>\n"
        "     <TD>%d</TD><TD>%d</TD><TD>%.6g</TD><TD>%d</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
        "    </TR>\n",
	e+1,
	nloaded_total, nloaded_min, nloaded_mean, nloaded_max,
	naccepted_total, naccepted_min, naccepted_mean, naccepted_max,
	minrad_min, minrad_mean, minrad_max,
	sampling_min, sampling_mean, sampling_max,
	chi2_min, chi2_mean, chi2_max,
	fwhm_min, fwhm_mean, fwhm_max,
	fwhm_wcs_min, fwhm_wcs_mean, fwhm_wcs_max,
	ellipticity_min, ellipticity_mean, ellipticity_max,
	ellipticity1_min, ellipticity1_mean, ellipticity1_max,
	ellipticity2_min, ellipticity2_mean, ellipticity2_max,
	beta_min, beta_mean, beta_max,
	residuals_min, residuals_mean, residuals_max,
	pffwhm_min, pffwhm_mean, pffwhm_max,
	pffwhm_wcs_min, pffwhm_wcs_mean, pffwhm_wcs_max,
	pfellipticity_min, pfellipticity_mean, pfellipticity_max,
	pfellipticity1_min, pfellipticity1_mean, pfellipticity1_max,
	pfellipticity2_min, pfellipticity2_mean, pfellipticity2_max,
	pfbeta_min, pfbeta_mean, pfbeta_max,
	pfresiduals_min, pfresiduals_mean, pfresiduals_max,
	symresiduals_min, symresiduals_mean, symresiduals_max,
	noiseqarea_min, noiseqarea_mean, noiseqarea_max,
	pixscale_wcs_min, pixscale_wcs_mean, pixscale_wcs_max);
    }

  fprintf(file, "   </TABLEDATA></DATA>\n");
  fprintf(file, "  </TABLE>\n");

/* Warnings */
  fprintf(file, "  <TABLE ID=\"Warnings\" name=\"Warnings\">\n");
  fprintf(file,
	"   <DESCRIPTION>%s warnings (limited to the last %d)</DESCRIPTION>\n",
	BANNER, WARNING_NMAX);
  fprintf(file, "   <FIELD name=\"Date\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta;time.event.end\"/>\n");
  fprintf(file, "   <FIELD name=\"Time\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta;time.event.end\"/>\n");
  fprintf(file, "   <FIELD name=\"Msg\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta\"/>\n");
  fprintf(file, "   <DATA><TABLEDATA>\n");
  for (str = warning_history(); *str; str = warning_history())
    fprintf(file, "    <TR><TD>%10.10s</TD><TD>%8.8s</TD><TD>%s</TD></TR>\n",
	str, str+11, str+22);
  fprintf(file, "   </TABLEDATA></DATA>\n");
  fprintf(file, "  </TABLE>\n");

/* Configuration file */
  fprintf(file, "  <RESOURCE ID=\"Config\" name=\"Config\">\n");
  fprintf(file, "   <DESCRIPTION>%s configuration</DESCRIPTION>\n", BANNER);
  fprintf(file,
	"   <PARAM name=\"Command_Line\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"obs.param\" value=\"%s",
	prefs.command_line[0]);
  for (n=1; n<prefs.ncommand_line; n++)
    fprintf(file, " %s", prefs.command_line[n]);
  fprintf(file, "\"/>\n");
  fprintf(file,
	"   <PARAM name=\"Prefs_Name\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"obs.param;meta.file\" value=\"%s\"/>\n",
	prefs.prefs_name);

  if (!error)
    {
/*-- PSF model */
    write_xmlconfigparam(file, "Basis_Type", "", "meta.code","%s");
    write_xmlconfigparam(file, "Basis_Number", "", "meta.number","%d");
    write_xmlconfigparam(file, "Basis_Name", "", "meta.id;meta.file","%s");
    write_xmlconfigparam(file, "Basis_Scale", "", "arith.factor","%.6g");
    write_xmlconfigparam(file, "NewBasis_Type", "", "meta.code","%s");
    write_xmlconfigparam(file, "NewBasis_Number", "", "meta.number","%d");
    write_xmlconfigparam(file, "PSF_Sampling", "pix",
		"arith.factor;instr.pixel;inst.det.psf","%.6g");
    write_xmlconfigparam(file, "PSF_Accuracy", "",
		"obs.param;phot.flux.sb;arith.ratio;instr.det.psf","%.6g");
    write_xmlconfigparam(file, "PSF_Size", "pix",
		"meta.number;instr.pixel;instr.det.psf","%d");
    write_xmlconfigparam(file, "Center_Keys", "",
		"meta.id;src;instr.det.psf", "%s");
    write_xmlconfigparam(file, "PSF_Recenter", "", "meta.code","%c");
    write_xmlconfigparam(file, "PhotFlux_Key", "",
		"meta.id;src;instr.det.psf", "%s");
    write_xmlconfigparam(file, "PhotFluxErr_Key", "",
		"meta.id;src;instr.det.psf", "%s");
    write_xmlconfigparam(file, "MEF_Type", "", "meta.code","%s");

/*-- PSF dependencies */
    write_xmlconfigparam(file, "PSFVar_Keys", "",
		"meta.id;src;instr.det.psf", "%s");
    write_xmlconfigparam(file, "PSFVar_Groups", "",
		"meta.id;stat.fit.param;instr.det.psf", "%d");
    write_xmlconfigparam(file, "PSFVar_Degrees", "",
		"stat.fit.param;instr.det.psf", "%d");
    write_xmlconfigparam(file, "PSFVar_NSnap", "",
		"stat.fit.param;instr.det.psf", "%d");
    write_xmlconfigparam(file, "HiddenMEF_Type", "", "meta.code","%s");
    write_xmlconfigparam(file, "Stability_Type", "", "meta.code","%s");

/*-- Sample selection */
    write_xmlconfigparam(file, "Sample_AutoSelect", "", "meta.code","%c");
    write_xmlconfigparam(file, "SampleVar_Type", "", "meta.code","%s");
    write_xmlconfigparam(file, "Sample_FWHMRange", "pix",
		"phys.size.diameter;instr.det.psf","%.6g");
    write_xmlconfigparam(file, "Sample_Variability", "",
		"instr.det.psf;arith.ratio","%.6g");
    write_xmlconfigparam(file, "Sample_MinSN", "", "stat.snr;stat.min","%.6g");
    write_xmlconfigparam(file, "Sample_MaxEllip", "",
		"src.ellipticity;stat.max","%.6g");
    write_xmlconfigparam(file, "Sample_FlagMask", "", "meta.code","%d");
    write_xmlconfigparam(file, "BadPixel_Filter", "", "meta.code","%c");
    write_xmlconfigparam(file, "BadPixel_NMax", "",
		"meta.number;instr.pixel;stat.max","%d");

/*-- Homogenisation kernel */
    write_xmlconfigparam(file, "HomoBasis_Type", "", "meta.code","%s");
    write_xmlconfigparam(file, "HomoBasis_Number", "", "meta.number","%d");
    write_xmlconfigparam(file, "HomoBasis_Scale", "", "arith.factor","%.6g");
    write_xmlconfigparam(file, "HomoPSF_Params", "", "stat.param","%.6g");
    write_xmlconfigparam(file, "HomoKernel_Dir", "", "meta.id;meta.file","%s");
    write_xmlconfigparam(file, "HomoKernel_Suffix", "",
		"meta.id;meta.file","%s");

/*-- Check-plots */
    write_xmlconfigparam(file, "CheckPlot_Dev", "", "meta.code", "%s");
    write_xmlconfigparam(file, "CheckPlot_Res", "", "meta.number;meta", "%d");
    write_xmlconfigparam(file, "CheckPlot_AntiAlias", "", "meta.code", "%c");
    write_xmlconfigparam(file, "CheckPlot_Type", "", "meta.code", "%s");
    write_xmlconfigparam(file, "CheckPlot_Name", "", "meta.id;meta.file", "%s");

/*-- Check Images --*/
    write_xmlconfigparam(file, "CheckImage_Type", "", "meta.code", "%s");
    write_xmlconfigparam(file, "CheckImage_Name", "",
		"meta.id;meta.file;meta.fits", "%s");
    write_xmlconfigparam(file, "CheckImage_Cube", "", "meta.code","%c");

/*-- Miscellaneous */
    write_xmlconfigparam(file, "PSF_Dir", "", "meta.id;meta.file","%s");
    write_xmlconfigparam(file, "PSF_Suffix", "", "meta.id;meta.file","%s");
    write_xmlconfigparam(file, "Verbose_Type", "", "meta.code","%s");
    write_xmlconfigparam(file, "Write_XML", "", "meta.code","%s");
    write_xmlconfigparam(file, "NThreads", "",
		"meta.number;meta.software", "%d");
    }

  fprintf(file, "  </RESOURCE>\n");
  fprintf(file, " </RESOURCE>\n");

  return RETURN_OK;
  }


/****** write_xmlerror ******************************************************
PROTO	int	write_xmlerror(char *error)
PURPOSE	Save meta-data to a simplified XML file in case of a catched error
INPUT	a character string.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	26/07/2008
 ***/
void	write_xmlerror(char *filename, char *error)
  {
   FILE			*file;
   int			pipe_flag;

  pipe_flag = 0;
  if (!strcmp(filename, "STDOUT"))
    {
    file = stdout;
    pipe_flag = 1;
    }
  else if (!(file = fopen(filename, "w")))
    return;

  write_xml_header(file);
  write_xml_meta(file, error);

  fprintf(file, "</RESOURCE>\n");
  fprintf(file, "</VOTABLE>\n");

  if (!pipe_flag)
    fclose(file);

  return;
  }


/****** write_xmlconfigparam **************************************************
PROTO	int write_xmlconfigparam(FILE *file, char *name, char *unit,
		char *ucd, char *format)
PURPOSE	Write to a VO-table the configuration parameters.
INPUT	Output stream (file) pointer,
	Name of the parameter keyword,
	unit,
	UCD string,
	printf() format to use in "value".
OUTPUT	RETURN_OK if the keyword exists, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	26/02/2014
 ***/
int	write_xmlconfigparam(FILE *file, char *name, char *unit,
		 char *ucd, char *format)
  {
   char		value[MAXCHAR], uunit[MAXCHAR];
   int		i,j,n;

  for (i=0; key[i].name[0] && cistrcmp(name, key[i].name, FIND_STRICT); i++);
  if (!key[i].name[0])
    return RETURN_ERROR;

  if (*unit)
    sprintf(uunit, " unit=\"%s\"", unit);
  else
    *uunit = '\0';
  switch(key[i].type)
    {
    case P_FLOAT:
      sprintf(value, format, *((double *)key[i].ptr));
      fprintf(file, "   <PARAM name=\"%s\"%s datatype=\"double\""
	" ucd=\"%s\" value=\"%s\"/>\n",
	name, uunit, ucd, value);
      break;
    case P_FLOATLIST:
      n = *(key[i].nlistptr);
      if (n)
        {
        sprintf(value, format, ((double *)key[i].ptr)[0]);
        fprintf(file, "   <PARAM name=\"%s\"%s datatype=\"double\""
		" arraysize=\"%d\" ucd=\"%s\" value=\"%s",
		name, uunit, n, ucd, value);
        for (j=1; j<n; j++)
          {
          sprintf(value, format, ((double *)key[i].ptr)[j]);
          fprintf(file, " %s", value);
          }
        fprintf(file, "\"/>\n");
        }
      else
        fprintf(file, "   <PARAM name=\"%s\"%s datatype=\"double\""
		" ucd=\"%s\" value=\"\"/>\n",
		name, uunit, ucd);
      break;
    case P_INT:
      sprintf(value, format, *((int *)key[i].ptr));
      fprintf(file, "   <PARAM name=\"%s\"%s datatype=\"int\""
	" ucd=\"%s\" value=\"%s\"/>\n",
	name, uunit, ucd, value);
      break;
    case P_INTLIST:
      n = *(key[i].nlistptr);
      if (n)
        {
        sprintf(value, format, ((int *)key[i].ptr)[0]);
        fprintf(file, "   <PARAM name=\"%s\"%s datatype=\"int\""
		" arraysize=\"%d\" ucd=\"%s\" value=\"%s",
		name, uunit, n, ucd, value);
        for (j=1; j<n; j++)
          {
          sprintf(value, format, ((int *)key[i].ptr)[j]);
          fprintf(file, " %s", value);
          }
        fprintf(file, "\"/>\n");
        }
      else
        fprintf(file, "   <PARAM name=\"%s\"%s datatype=\"double\""
		" ucd=\"%s\" value=\"\"/>\n",
		name, uunit, ucd);
      break;
    case P_BOOL:
      sprintf(value, "%c", *((int *)key[i].ptr)? 'T':'F');
      fprintf(file, "   <PARAM name=\"%s\" datatype=\"boolean\""
	" ucd=\"%s\" value=\"%s\"/>\n",
	name, ucd, value);
      break;
    case P_BOOLLIST:
      n = *(key[i].nlistptr);
      if (n)
        {
        sprintf(value, "%c", ((int *)key[i].ptr)[0]? 'T':'F');
        fprintf(file, "   <PARAM name=\"%s\" datatype=\"boolean\""
		" arraysize=\"%d\" ucd=\"%s\" value=\"%s",
		name, n, ucd, value);
        for (j=1; j<n; j++)
          {
          sprintf(value, "%c", ((int *)key[i].ptr)[j]? 'T':'F');
          fprintf(file, " %s", value);
          }
        fprintf(file, "\"/>\n");
        }
      else
        fprintf(file, "   <PARAM name=\"%s\" datatype=\"boolean\""
		" ucd=\"%s\" value=\"\"/>\n",
		name, ucd);
      break;
    case P_STRING:
      sprintf(value, "%s", (char *)key[i].ptr);
      fprintf(file, "   <PARAM name=\"%s\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"%s\" value=\"%s\"/>\n",
	name, ucd, *value? value: " ");
      break;
    case P_STRINGLIST:
      n = *(key[i].nlistptr);
      if (n)
        {
        sprintf(value, "%s", ((char **)key[i].ptr)[0]);
        fprintf(file, "   <PARAM name=\"%s\" datatype=\"char\""
		" arraysize=\"*\" ucd=\"%s\" value=\"%s",
		name, ucd, *value? value: " ");
        for (j=1; j<n; j++)
          {
          sprintf(value, "%s", ((char **)key[i].ptr)[j]);
          fprintf(file, ",%s", *value? value: " ");
          }
        fprintf(file, "\"/>\n");
        }
      else
        fprintf(file, "   <PARAM name=\"%s\" datatype=\"char\""
		" arraysize=\"*\" ucd=\"%s\" value=\"\"/>\n",
		name, ucd);
      break;
    case P_KEY:
      sprintf(value, "%s", key[i].keylist[*((int *)key[i].ptr)]);
      fprintf(file, "   <PARAM name=\"%s\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"%s\" value=\"%s\"/>\n",
	name, ucd, value);
      break;
    case P_KEYLIST:
      n = *(key[i].nlistptr);
      if (n)
        {
        sprintf(value, "%s", key[i].keylist[((int *)key[i].ptr)[0]]);
        fprintf(file, "   <PARAM name=\"%s\" datatype=\"char\""
		" arraysize=\"*\" ucd=\"%s\" value=\"%s",
		name, ucd, value);
        for (j=1; j<n; j++)
          {
          sprintf(value, "%s", key[i].keylist[((int *)key[i].ptr)[j]]);
          fprintf(file, ",%s", value);
          }
        fprintf(file, "\"/>\n");
        }
      else
        fprintf(file, "   <PARAM name=\"%s\" datatype=\"char\""
		" arraysize=\"*\" ucd=\"%s\" value=\"\"/>\n",
		name, ucd);
      break;
    default:
        error(EXIT_FAILURE, "*Internal Error*: Type Unknown",
		" in write_xmlconfigparam()");
    }

  return RETURN_OK;
  }

