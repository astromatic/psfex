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
*	Last modify:	01/03/2007
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
EXECUTABLE " catalog [catalog2 ...][-c <config_file>][-<keyword> <value>]\n" \
"> to dump a default configuration file: " EXECUTABLE " -d \n" \
"> to dump a default extended configuration file: " EXECUTABLE " -dd \n"

/********************************** main ************************************/

int main(int argc, char *argv[])
  {
   int		a,ab,na, narg, opt, opt2;
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
  prefs.command_line = argv;
  prefs.ncommand_line = argc;
  ab=na=0;
  narg = 0;
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

  NFPRINTF(OUTPUT, "");
  NPRINTF(OUTPUT, "> All done (in %.0f s)\n", prefs.time_diff);

  exit(EXIT_SUCCESS);
  return 0;
  }
