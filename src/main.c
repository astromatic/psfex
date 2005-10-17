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
*	Last modify:	31/10/2003
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<ctype.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"types.h"
#include	"globals.h"
#include	"fits/fitscat.h"
#include	"prefs.h"

#define		SYNTAX \
EXECUTABLE " catalog [catalog2 ...][-c <config_file>]" \
	"[-<keyword> <value>]\n" \
	"> or, to dump a default configuration file:\n" \
	"> " EXECUTABLE " -d \n"

/********************************** main ************************************/

int main(int argc, char *argv[])
  {
   int		a,ab,na, narg, opt;
   char		**argkey, **argval;

  if (argc<2)
    {
    fprintf(OUTPUT, "\n         %s  Version %s (%s)\n", BANNER,MYVERSION,DATE);
    fprintf(OUTPUT, "\nby %s\n", COPYRIGHT);
    fprintf(OUTPUT, "visit %s\n", WEBSITE);
    error(EXIT_SUCCESS, "SYNTAX: ", SYNTAX);
    }

  QMALLOC(argkey, char *, argc);
  QMALLOC(argval, char *, argc);

/* Default parameters */
  ab=na=0;
  narg = 0;
  strcpy(prefs.prefs_name, "default.psfex");

  for (a=1; a<argc; a++)
    {
    if (*(argv[a]) == '-')
      {
      opt = (int)argv[a][1];
      if (strlen(argv[a])<3 || opt == '-')
        {
        if (opt == '-')
          opt = (int)tolower((int)argv[a][2]);
        switch(opt)
          {
          case 'c':
            if (a<(argc-1))
              strcpy(prefs.prefs_name, argv[++a]);
            break;
          case 'd':
            dumpprefs();
            exit(EXIT_SUCCESS);
            break;
          case 'v':
            printf("%s version %s (%s)\n", BANNER,MYVERSION,DATE);
            exit(0);
            break;
          case 'h':
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
      for(ab = a; (a<argc) && (*argv[a]!='-'); a++);
      na = (a--) - ab;
      }
    }

  if (!na)
    error(EXIT_SUCCESS,"SYNTAX: ", SYNTAX);

  readprefs(prefs.prefs_name, argkey, argval, narg);
  useprefs();
  free(argkey);
  free(argval);


  makeit(argv+ab, na);

  NFPRINTF(OUTPUT, "All done");
  NPRINTF(OUTPUT, "\n");

  exit(EXIT_SUCCESS);
  return 0;
  }
