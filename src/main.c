 /*
 				main.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	PSFEx
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	parsing and main loop.
*
*	Last modify:	19/06/2009
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifdef HAVE_PLPLOT
#include PLPLOT_H
#endif

#include "define.h"
#include "types.h"
#include "globals.h"
#include "fits/fitscat.h"
#include "prefs.h"
#include "cplot.h"

#define		SYNTAX \
EXECUTABLE " catalog1 [catalog2,...][@catalog_list1 [@catalog_list2 ...]]\n" \
"\t\t[-c <config_file>][-<keyword> <value>]\n" \
"> to dump a default configuration file: " EXECUTABLE " -d \n" \
"> to dump a default extended configuration file: " EXECUTABLE " -dd \n"

extern void	(*myplparseopts)(int *p_argc, const char **argv, PLINT mode);

extern const char       notokstr[];

/********************************** main ************************************/

int main(int argc, char *argv[])
  {
   char		**argkey, **argval,
		*str,*listbuf;
   int		a, narg, nim, ntok, opt, opt2;

#ifdef HAVE_SETLINEBUF
/* flush output buffer at each line */
  setlinebuf(stderr);
#endif

  if (argc<2)
    {
    fprintf(OUTPUT, "\n         %s  Version %s (%s)\n", BANNER,MYVERSION,DATE);
    fprintf(OUTPUT, "\nby %s\n", COPYRIGHT);
    fprintf(OUTPUT, "visit %s\n", WEBSITE);
    error(EXIT_SUCCESS, "SYNTAX: ", SYNTAX);
    }

#ifdef HAVE_PLPLOT
  cplot_fixplplot();
  if (argc>2)
    myplparseopts(&argc, (const char **)argv, PL_PARSE_SKIP);
#endif

  QMALLOC(argkey, char *, argc);
  QMALLOC(argval, char *, argc);

/* Default parameters */
  prefs.command_line = argv;
  prefs.ncommand_line = argc;
  narg = nim = 0;
  listbuf = (char *)NULL;
  strcpy(prefs.prefs_name, "default.psfex");

  for (a=1; a<argc; a++)
    {
    if (*(argv[a]) == '-')
      {
      opt = (int)argv[a][1];
      if (strlen(argv[a])<4 || opt == '-')
        {
        opt2 = (int)tolower((int)argv[a][2]);
        if (opt == '-')
          {
          opt = opt2;
          opt2 = (int)tolower((int)argv[a][3]);
          }
        switch(opt)
          {
          case 'c':
            if (a<(argc-1))
              strcpy(prefs.prefs_name, argv[++a]);
            break;
          case 'd':
            dumpprefs(opt2=='d' ? 1 : 0);
            exit(EXIT_SUCCESS);
            break;
          case 'v':
            printf("%s version %s (%s)\n", BANNER,MYVERSION,DATE);
            exit(EXIT_SUCCESS);
            break;
          case 'h':
            fprintf(OUTPUT, "\nSYNTAX: %s", SYNTAX);
#ifdef HAVE_PLPLOT
            fprintf(OUTPUT, "\nPLPLOT-specific options:\n");
            myplparseopts(&argc, (const char **)argv, PL_PARSE_SKIP);
#endif
            exit(EXIT_SUCCESS);
            break;
          default:
            error(EXIT_SUCCESS,"SYNTAX: ", SYNTAX);
          }
        }
      else
        {
/*------ Config parameters */
        argkey[narg] = &argv[a][1];
        argval[narg++] = argv[++a];
        }       
      }
    else
      {
/*---- The input image filename(s) */
      for(; (a<argc) && (*argv[a]!='-'); a++)
        {
        str = (*argv[a] == '@'? listbuf=list_to_str(argv[a]+1) : argv[a]);
        for (ntok=0; (str=strtok(ntok?NULL:str, notokstr)); nim++,ntok++)
          if (nim<MAXFILE)
            prefs.incat_name[nim] = str;
          else
            error(EXIT_FAILURE, "*Error*: Too many input catalogues: ", str);
        }
      a--;
      }
    }

  prefs.ncat = nim;

  readprefs(prefs.prefs_name, argkey, argval, narg);
  useprefs();
  free(argkey);
  free(argval);

  makeit();

  free(listbuf);

  NFPRINTF(OUTPUT, "");
  NPRINTF(OUTPUT, "> All done (in %.0f s)\n", prefs.time_diff);

  exit(EXIT_SUCCESS);
  }
