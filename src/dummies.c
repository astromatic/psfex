#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "define.h"
#include "sample.h"

struct cat; typedef struct cat catstruct;
struct field; typedef struct field fieldstruct;
struct key; typedef struct key keystruct;
struct tab; typedef struct tab tabstruct;
struct wcs; typedef struct wcs wcsstruct;

typedef enum {H_TYPE} h_type;
typedef enum {T_TYPE} t_type;
typedef enum {ACCESS_TYPE_T} access_type_t;
typedef float PIXTYPE;

const int t_size[] = {0};
time_t thetime, thetime2;

int
add_key(keystruct *key, tabstruct *tab, int pos)
{
   abort();
   return -1;
}

int
addkeywordto_head(tabstruct *tab, char *keyword, char *comment)
{
   abort();
   return -1;
}

int
blank_keys(tabstruct *tab)
{
   abort();
   return -1;
}

void
end_wcs(wcsstruct *wcs)
{
   if (wcs != NULL) {
      abort();
   }
}

int
fitsread(char *fitsbuf, char *keyword, void *ptr, h_type htype,
	 t_type ttype)
{
   abort();
   return -1;
}

int
fitswrite(char *fitsbuf, char *keyword, void *ptr, h_type htype,
	  t_type ttype)
{
   abort();
   return -1;
}

int
init_cat(catstruct *cat)
{
   abort();
   return -1;
}

catstruct *
read_cat(char *filename)
{
   abort();
   return NULL;
}

keystruct *
read_key(tabstruct *tab, char *keyname)
{
   abort();
   return NULL;
}

keystruct *
new_key(char *keyname)
{
   abort();
   return NULL;
}

tabstruct *
new_tab(char *tabname)
{
   abort();
   return NULL;
}

wcsstruct *
read_wcs(tabstruct *tab)
{
   abort();
   return NULL;
}

double
wcs_dist(wcsstruct *wcs, double *wcspos1, double *wcspos2)
{
   abort();
   return -1.0;
}

double
wcs_scale(wcsstruct *wcs, double *pixpos)
{
#if defined(NAN)
   return NAN;
#else
   return -1.0;
#endif
}
   
void
readbasic_head(tabstruct *tab)
{
   abort();
}

catstruct *
new_cat(int ncat)
{
   abort();
   return NULL;
}

void
free_cat(catstruct **cat, int ncat)
{
   abort();
}

void
free_tab(tabstruct *tab)
{
   abort();
}

int
open_cat(catstruct *cat, access_type_t at)
{
   abort();
   return -1;
}

int
prim_head(tabstruct *tab)
{
   abort();
   return -1;
}

void
read_body(tabstruct *tab, PIXTYPE *ptr, size_t size)
{
   abort();
}

void
save_tab(catstruct *cat, tabstruct *tab)
{
   abort();
}

/*****************************************************************************/
/*
 * These are not dummies.  They replace psfex routines
 */
typedef enum {CHECKENUM} checkenum;

void
check_write(fieldstruct *field, setstruct *set,
	    char *checkname, checkenum checktype,
	    int ext, int next, int cubeflag)
{
    ;
}

setstruct *
load_samples(char **filenames, int catindex, int ncat, int ext,
             int next, contextstruct *context)
{
    /*
     * The C version of this is called two ways:
     *   catindex == 0, ncat == ncat            Read all catalogues
     *   catindex == c, ncat == 1               Read only catalogue c
     */
   setstruct *completeSet = (setstruct *)(filenames[catindex + 0]);
    /*
     * Make a new set, which may be a subset of the completeSet
     */
    setstruct *set = init_set(context);
    set->fwhm = completeSet->fwhm;
    for (int i = 0; i != completeSet->vigdim; ++i) {
        set->vigsize[i] = completeSet->vigsize[i];
    }
    for (int i = 0; i != completeSet->ncontext; ++i) {
        strcpy(set->contextname[i], completeSet->contextname[i]);
        set->contextoffset[i] = completeSet->contextoffset[i];
        set->contextscale[i] = completeSet->contextscale[i];
    }
    /*
     * Count how many samples we'll be including
     */
    int nsample_keep = 0;
    for (int i = 0; i != ncat; ++i) {
        setstruct *s = (setstruct *)(filenames[catindex + i]);
        for (int j = 0; j != completeSet->nsample; ++j) {
            samplestruct const *samp = s->sample[j];
            if (ext == ALL_EXTENSIONS || ext == samp->extindex) {
                ++nsample_keep;
            }
        }
    }

    set->samples_owner = 0;
    malloc_samples(set, nsample_keep);
    for (int i = 0; i != ncat; ++i) {
       setstruct *s = (setstruct *)(filenames[catindex + i]);
        for (int j = 0; j != completeSet->nsample; ++j) {
            samplestruct *samp = s->sample[j];
            if (ext == ALL_EXTENSIONS || ext == samp->extindex) {
                set->sample[set->nsample++] = samp;
            }
        }
    }

    return set;
}
