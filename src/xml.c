 /*
				xml.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	XML logging.
*
*	Last modify:	28/03/2007
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

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
#include "key.h"
#include "prefs.h"
#include "psf.h"
#include "xml.h"

extern time_t		thetime,thetime2;	/* from makeit.c */
extern pkeystruct	key[];			/* from preflist.h */
extern char		keylist[][32];		/* from preflist.h */
 
psfstruct		**psf_xml;
int			*nfield_xml;
int			nxml, nxmlmax;


/****** init_xml ************************************************************
PROTO	int init_xml(int next)
PURPOSE	Initialize a set of meta-data kept in memory before being written to the
	XML file
INPUT	Number of extensions.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	02/03/2007
 ***/
int	init_xml(int next)
  {
  QMALLOC(psf_xml, psfstruct *, next);
  QMALLOC(nfield_xml, int, next);
  nxml = 0;
  nxmlmax = next;

  return EXIT_SUCCESS;
  }


/****** end_xml ************************************************************
PROTO	void end_xml(void)
PURPOSE	Free the set of meta-data kept in memory.
INPUT	-.
OUTPUT	.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	02/03/2007
 ***/
void	end_xml(void)
  {
  free(psf_xml);
  free(nfield_xml);

  return;
  }


/****** update_xml ***********************************************************
PROTO	int update_xml(psfstruct *psf, int nfields)
PURPOSE	Update a set of meta-data kept in memory before being written to the
	XML file
INPUT	Pointer to the current PSF,
	number of fields loaded.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	Global preferences are used.
AUTHOR	E. Bertin (IAP)
VERSION	02/03/2007
 ***/
int	update_xml(psfstruct *psf, int nfield)
  {
  psf_xml[nxml] = psf;
  nfield_xml[nxml] = nfield;
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
VERSION	06/10/2006
 ***/
int	write_xml_header(FILE *file)
  {
   char		sysname[16];

  fprintf(file, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
  fprintf(file, "<?xml-stylesheet type=\"text/xsl\" href=\"%s\"?>\n",
	prefs.xsl_name);
  fprintf(file, "<VOTABLE "
	"xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
	"xsi:noNamespaceSchemaLocation="
	"\"http://www.ivoa.net/xml/VOTable/v1.1\">\n");
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
VERSION	28/03/2007
 ***/
int	write_xml_meta(FILE *file, char *error)
  {
   psfstruct		*psf;
   struct tm		*tm;
   char			*pspath,*psuser, *pshost, *str;
   double		fwhm_best,fwhm_worse, elongation_best,elongation_worse,
			beta_best,beta_worse, residuals_best,residuals_worse,
			temp;
   int			i,d,n, nmed, ntot;

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

/* Meta-data for the PSF */
  fprintf(file, "  <TABLE ID=\"PSF\" name=\"PSF\">\n");
  fprintf(file, "   <DESCRIPTION>Metadata and stats about the PSF gathered by "
	"%s</DESCRIPTION>\n", BANNER);
  fprintf(file, "   <!-- Nextensions may be 0"
	" if an error occurred early in the processing -->\n");
  fprintf(file, "   <PARAM name=\"NExtensions\" datatype=\"int\""
	" ucd=\"meta.number;meta.dataset\" value=\"%d\"/>\n", nxmlmax);
  fprintf(file, "   <!-- CurrExtension may differ from Nextensions"
	" if an error occurred -->\n");
  fprintf(file, "   <PARAM name=\"CurrExtension\" datatype=\"int\""
	" ucd=\"meta.number;meta.dataset\" value=\"%d\"/>\n", nxml);
  fprintf(file, "   <FIELD name=\"NStars_Loaded\" datatype=\"int\""
	" ucd=\"meta.number;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"NStars_Accepted\" datatype=\"int\""
	" ucd=\"meta.number;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_FromFluxRadius\" unit=\"pix\""
	" datatype=\"float\" ucd=\"phys.size.diameter;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Sampling\" unit=\"pix\" datatype=\"float\""
	" ucd=\"arith.factor;instr.pixel;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Chi2\" datatype=\"float\""
	" ucd=\"stat.fit.chi2;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"NSnapshots\" datatype=\"int\""
	" arraysize=\"%d\" ucd=\"meta.number;meta.dataset\"/>\n",
	nxml? psf_xml[0]->poly->ndim: 1);
  fprintf(file, "   <FIELD name=\"FWHM\" unit=\"pix\""
	" datatype=\"float\" ucd=\"phys.size.diameter;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_Best\" unit=\"pix\" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"FWHM_Worse\" unit=\"pix\" datatype=\"float\""
	" ucd=\"phys.size.diameter;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Elongation\" datatype=\"float\""
	" ucd=\"phys.size.axisRatio;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Elongation_Best\" datatype=\"float\""
	" ucd=\"phys.size.axisRatio;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Elongation_Worse\" datatype=\"float\""
	" ucd=\"phys.size.axisRatio;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"MoffatBeta\" datatype=\"float\""
	" ucd=\"stat.param;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"MoffatBeta_Best\" datatype=\"float\""
	" ucd=\"stat.param;stat.max;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"MoffatBeta_Worse\" datatype=\"float\""
	" ucd=\"stat.param;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Residuals\" datatype=\"float\""
	" ucd=\"stat.fit.residual;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Residuals_Best\" datatype=\"float\""
	" ucd=\"stat.fit.residual;stat.min;instr.det.psf\"/>\n");
  fprintf(file, "   <FIELD name=\"Residuals_Worse\" datatype=\"float\""
	" ucd=\"stat.fit.residual;stat.max;instr.det.psf\"/>\n");

  fprintf(file, "   <DATA><TABLEDATA>\n");
  for (n=0; n<nxml; n++)
    {
    psf = psf_xml[n];
    fprintf(file, "    <TR>\n"
        "     <TD>%d</TD><TD>%d</TD><TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%d",
	psf->samples_loaded,
	psf->samples_accepted,
	psf->fwhm,
        psf->pixstep,
	psf->chi2,
	psf->poly->ndim? prefs.context_nsnap : 1);
    ntot =  1;
    nmed = 0;
    for (d=0; d<psf_xml[0]->poly->ndim; d++)
      {
      fprintf(file, " %d", prefs.context_nsnap);
      ntot *= prefs.context_nsnap;
      nmed += nmed*prefs.context_nsnap + (prefs.context_nsnap-1)/2;
      }
    fwhm_best = elongation_best = beta_worse = residuals_best = BIG;
    fwhm_worse = elongation_worse = beta_best = residuals_worse = -BIG;
    for (i=0; i<ntot; i++)
      {
      if ((temp=sqrt(psf->moffat[i].fwhm_min*psf->moffat[i].fwhm_max))
	< fwhm_best)
        fwhm_best = temp;
      if (temp > fwhm_worse)
        fwhm_worse = temp;
      if ((temp=psf->moffat[i].fwhm_max / psf->moffat[i].fwhm_min)
	< elongation_best)
        elongation_best = temp;
      if (temp > elongation_worse)
        elongation_worse = temp;
      if (psf->moffat[i].beta > beta_best)
        beta_best = psf->moffat[i].beta;
      if (psf->moffat[i].beta < beta_worse)
        beta_worse = psf->moffat[i].beta;
      if (psf->moffat[i].residuals < residuals_best)
        residuals_best = psf->moffat[i].residuals;
      if (psf->moffat[i].residuals > residuals_worse)
        residuals_worse = psf->moffat[i].residuals;
      }

    fprintf(file, "</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
	"     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
        "    </TR>\n",
	sqrt(psf->moffat[nmed].fwhm_min*psf->moffat[nmed].fwhm_max),
	fwhm_best,
	fwhm_worse,
	psf->moffat[nmed].fwhm_max/psf->moffat[nmed].fwhm_min,
	elongation_best,
	elongation_worse,
	psf->moffat[nmed].beta,
	beta_best,
	beta_worse,
	psf->moffat[nmed].residuals,
	residuals_best,
	residuals_worse);
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
    write_xmlconfigparam(file, "PSF_Name", "",
		"meta.id;meta.file","%s");
    write_xmlconfigparam(file, "PSF_Accuracy", "",
		"obs.param;phot.flux.sb;arith.ratio","%.6g");
    write_xmlconfigparam(file, "PSF_NSuper", "pix",
		"meta.number;instr.pixel","%d");
    write_xmlconfigparam(file, "PSF_Sampling", "pix",
		"arith.factor;instr.pixel","%.6g");
    write_xmlconfigparam(file, "PSF_Size", "pix",
		"meta.number;instr.pixel","%d");
    write_xmlconfigparam(file, "PSF_Recenter", "",
		"meta.code","%c");

/*-- Sample selection */
    write_xmlconfigparam(file, "PSF_AutoSelect", "",
		"meta.code","%c");
    write_xmlconfigparam(file, "PSF_FWHMRange", "pix",
		"phys.size.diameter;instr.det.psf","%.6g");
    write_xmlconfigparam(file, "PSF_Variability", "",
		"instr.det.psf;arith.ratio","%.6g");
    write_xmlconfigparam(file, "PSF_MinSN", "",
		"stat.snr;stat.min","%.6g");
    write_xmlconfigparam(file, "PSF_MaxElong", "",
		"src.ellipticity;stat.max","%.6g");
    write_xmlconfigparam(file, "BadPixel_Filter", "",
		"meta.code","%c");
    write_xmlconfigparam(file, "BadPixel_NMax", "",
		"meta.number;instr.pixel;stat.max","%d");

/*-- PSF dependencies */
    write_xmlconfigparam(file, "Context_Keys", "",
		"meta.id;src", "%s");
    write_xmlconfigparam(file, "Context_Groups", "",
		"meta.id;stat.fit.param", "%d");
    write_xmlconfigparam(file, "Context_Degrees", "",
		"meta.id;stat.fit.param", "%d");
    write_xmlconfigparam(file, "CheckImage_Type", "",
		"meta.code", "%s");
    write_xmlconfigparam(file, "CheckImage_Name", "",
		"meta.id;meta.file;meta.fits", "%s");

/*-- Miscellaneous */
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
VERSION	26/07/2006
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
VERSION	06/10/2006
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
      sprintf(value, (char *)key[i].ptr);
      fprintf(file, "   <PARAM name=\"%s\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"%s\" value=\"%s\"/>\n",
	name, ucd, value);
      break;
    case P_STRINGLIST:
      n = *(key[i].nlistptr);
      if (n)
        {
        sprintf(value, ((char **)key[i].ptr)[0]);
        fprintf(file, "   <PARAM name=\"%s\" datatype=\"char\""
		" arraysize=\"*\" ucd=\"%s\" value=\"%s",
		name, ucd, value);
        for (j=1; j<n; j++)
          {
          sprintf(value, ((char **)key[i].ptr)[j]);
          fprintf(file, ",%s", value);
          }
        fprintf(file, "\"/>\n");
        }
      else
        fprintf(file, "   <PARAM name=\"%s\" datatype=\"char\""
		" arraysize=\"*\" ucd=\"%s\" value=\"\"/>\n",
		name, ucd);
      break;
    case P_KEY:
      sprintf(value, key[i].keylist[*((int *)key[i].ptr)]);
      fprintf(file, "   <PARAM name=\"%s\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"%s\" value=\"%s\"/>\n",
	name, ucd, value);
      break;
    case P_KEYLIST:
      n = *(key[i].nlistptr);
      if (n)
        {
        sprintf(value, key[i].keylist[((int *)key[i].ptr)[0]]);
        fprintf(file, "   <PARAM name=\"%s\" datatype=\"char\""
		" arraysize=\"*\" ucd=\"%s\" value=\"%s",
		name, ucd, value);
        for (j=1; j<n; j++)
          {
          sprintf(value, key[i].keylist[((int *)key[i].ptr)[j]]);
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

