#include <stdlib.h>
#include <stdio.h>
#include "netcdf.h"
#include "particleIO.h"

////////////////////////////////////////////////////////////////////////////////
// readdata: Read particle constants-of-motion data from a NetCDF file.
////////////////////////////////////////////////////////////////////////////////
void readdata(const char *fname, vvect **p, int *np)
{
  double *buffer;
  signed char   *sgnv;
  size_t len0, len1, start[2], count[2];
  int    dimids[2];
  int    ierr, ncid, ppid, muid, Eid, wtid, svid, ndims, ipart;

  /* Open the NetCDF file */
  if ((ierr = nc_open(fname, NC_NOWRITE, &ncid)) != NC_NOERR) {
    fprintf(stderr, "Error opening particle file %s: %s.\n",
	    fname, nc_strerror(ierr));
    exit(1);
  }

  /* Check for required variables */
  if ((ierr = nc_inq_varid(ncid, "pphi", &ppid)) != NC_NOERR) {
    fprintf(stderr, "pphi: %s\n", nc_strerror(ierr));
    exit(1);
  }
  if ((ierr = nc_inq_varid(ncid, "mu", &muid)) != NC_NOERR) {
    fprintf(stderr, "mu: %s\n", nc_strerror(ierr));
    exit(1);
  }
  if ((ierr = nc_inq_varid(ncid, "E", &Eid)) != NC_NOERR) {
    fprintf(stderr, "E: %s\n", nc_strerror(ierr));
    exit(1);
  }
  if ((ierr = nc_inq_varid(ncid, "weight", &wtid)) != NC_NOERR) {
    fprintf(stderr, "weight: %s\n", nc_strerror(ierr));
    exit(1);
  }
  if ((ierr = nc_inq_varid(ncid, "sgn_v", &svid)) != NC_NOERR) {
    fprintf(stderr, "sgn_v: %s\n", nc_strerror(ierr));
    exit(1);
  }

  /* Get variable dimensions */
  ierr = nc_inq_varndims(ncid, ppid, &ndims);
  if ((ndims != 2) || (ierr != NC_NOERR)) {
    fputs("Wrong variable dimension count.\n", stderr);
    exit(1);
  }
  ierr = nc_inq_vardimid(ncid, ppid, dimids);
  ierr = nc_inq_dimlen(ncid, dimids[0], &len0);
  ierr = nc_inq_dimlen(ncid, dimids[1], &len1);
  printf("%s contains %ld species, %ld particles\n", fname,
	 (long)len1, (long)len0);
  *np = (int)len0;

  /* Allocate storage */
  if ((*p = (vvect *)malloc(*np * sizeof(vvect))) == NULL) {
    fputs("Out of memory in readdata().\n", stderr);
    exit(10);
  }
  if ((buffer = (double *)malloc(*np * sizeof(double))) == NULL) {
    fputs("Out of memory in readdata().\n", stderr);
    exit(10);
  }
  fprintf(stderr, "Reading %e particles from %s...\n", (float)*np, fname);

  /* Read data */
  start[0] = start[1] = 0L;
  count[0] = len0;  count[1] = 1L;
  if ((ierr = nc_get_vara_double(ncid, ppid, start, count, buffer))
      != NC_NOERR) {
    fprintf(stderr, "Error reading pphi: %s\n", nc_strerror(ierr));
    exit(11);
  } else fputs("pphi read.\n", stderr);
  for (ipart=0; ipart<*np; ipart++) p[0][ipart].pphi = buffer[ipart];
  if ((ierr = nc_get_vara_double(ncid, muid, start, count, buffer))
      != NC_NOERR) {
    fprintf(stderr, "Error reading mu: %s\n", nc_strerror(ierr));
    exit(11);
  } else fputs("mu read.\n", stderr);
  for (ipart=0; ipart<*np; ipart++) p[0][ipart].mu = buffer[ipart];
  if ((ierr = nc_get_vara_double(ncid, Eid, start, count, buffer))
      != NC_NOERR) {
    fprintf(stderr, "Error reading E: %s\n", nc_strerror(ierr));
    exit(11);
  } else fputs("ke read.\n", stderr);
  for (ipart=0; ipart<*np; ipart++) p[0][ipart].ke = buffer[ipart];
  if ((ierr = nc_get_vara_double(ncid, wtid, start, count, buffer))
      != NC_NOERR) {
    fprintf(stderr, "Error reading weight: %s\n", nc_strerror(ierr));
    exit(11);
  } else fputs("weights read.\n", stderr);
  for (ipart=0; ipart<*np; ipart++) p[0][ipart].wt = buffer[ipart];
  free(buffer);

  sgnv = (signed char *)malloc(*np);
  if ((ierr = nc_get_vara_schar(ncid, svid, start, count, sgnv))
      != NC_NOERR) {
    fprintf(stderr, "Error reading sgn_v: %s\n", nc_strerror(ierr));
    exit(11);
  } else fputs("sgn_v read.\n", stderr);
  for (ipart=0; ipart<*np; ipart++) p[0][ipart].sig = (int)sgnv[ipart];
  free(sgnv);

  ierr = nc_close(ncid);
}

////////////////////////////////////////////////////////////////////////////////
// Comparison function for two particles, for quicksort.
////////////////////////////////////////////////////////////////////////////////
int mcompar(const void *f1, const void *f2)
{
  double mu1, mu2;
  int    s1 = ((vvect *)f1)->sig, s2 = ((vvect *)f2)->sig;

  if (s1==s2) {
    mu1 = ((vvect *)f1)->mu;  mu2 = ((vvect *)f2)->mu;
    //if (mu1 < mu2) return -1;
    //if (mu1 > mu2) return 1;
    //return 0;
    return (mu1 > mu2) - (mu1 < mu2);
  }

  return s1-s2;
}
