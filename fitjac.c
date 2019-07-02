#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "netcdf.h"
#include "spline_interface.h"

#define DEBUG
#define PPBIN 40

/* Type definitions */
typedef struct {
  double pphi, mu, ke, wt;
  int    sig;
} vvect;

/* Function prototypes */
void readdata(const char *, vvect **, int *);
int mcompar(const void *, const void *);
void fillbins(double ***arr, int nk, double cmin, double ew,
	      int np, double pmin, double pw,
	      vvect *part, int ji, int jf);
int killorphans(int **, const int, const int);
int threshold(int **, const int, const int, const int);
void gaussblur(double **, const int, const int, const double, const double);
int getlbound(vvect *, const int);

/******************************************************************************/
int main(int argc, char *argv[])
{
  FILE   *fp;
#ifdef DEBUG
  FILE   *fp2;
  char   f2name[64];
  double dtally;
#endif
  char outname[32], inname[96];
  vvect  *particle, *jacobian;
  spline **pspline;
  spline *kspline;
  double *ke, *pphi, *count;
  double *mubounds;
  double **binarr, **jacarr;
  const double sigma_ptcl = 1.334;
  const double sigma_jac = 1.75;
  const int minbins = 7;
  double  delta, ratio, pmin, pmax, cmin, emax;
  double  pbinwidth, ebinwidth, elvolinv;
  int     nparts, njac, nmubins, i, j, k, nkbins, npbins, npcoefs, nkcoefs;
  int     binstart, binstop, jbinstart, jbinstop, ko, kocount=0, lbound, jlbound;
  int     lmin=0, lmax, jlmin=0, jlmax, mu_eqpart=1;

  if (argc==1) {
    fputs("Enter file root.\n", stderr);
    exit(0);
  }

  /* Read all particle data from file */
  sprintf(inname, "%s_condensed.cdf", argv[1]);
  readdata(inname, &particle, &nparts);
  fprintf(stderr, "%d particles read; last = %le\t%le\t%le\t%d\n", nparts,
	  particle[nparts-1].pphi, particle[nparts-1].mu, particle[nparts-1].ke,
	  particle[nparts-1].sig);

  /* Read all jacobian data from file */
  sprintf(inname, "%s_jacobian.cdf", argv[1]);
  readdata(inname, &jacobian, &njac);
  fprintf(stderr, "%d particles read; last = %le\t%le\t%le\t%d\n", njac,
	  jacobian[nparts-1].pphi, jacobian[nparts-1].mu, jacobian[nparts-1].ke,
	  jacobian[nparts-1].sig);

  /* Calculate cmin (minimum ratio of ke to mu), emax */
  cmin = particle[0].ke / particle[0].mu;  emax = particle[0].ke;
  for (i=1; i<nparts; i++) {
    if (particle[i].ke > emax) emax = particle[i].ke;
    ratio = particle[i].ke / particle[i].mu;
    if (ratio < cmin) cmin = ratio;
  }
  fprintf(stderr, "cmin = %lf; emax = %le\n", cmin, emax);

  /* Sort the particles by sign of v||, increasing magnetic moment mu */
  fputs("Sorting by lambda, magnetic moment...\n", stderr);
  qsort(particle, nparts, sizeof(vvect), mcompar);
  lmax = lbound = getlbound(particle, nparts);
  fprintf(stderr, "%d have l<0 (%.1lf %%).\n", lbound, 100.0*lbound/nparts);
  fputs("Sorting Jacobian...\n", stderr);
  qsort(jacobian, njac, sizeof(vvect), mcompar);
  jlmax = jlbound = getlbound(jacobian, njac);
  fprintf(stderr, "%d have l<0 (%.1lf %%).\n", jlbound, 100.0*jlbound/njac);

  do { /* Once for lambda<0, once for lambda>0 */
    fprintf(stderr, " s=%d: mu ranges from %le to %le.\n", particle[lmin].sig,
	    particle[lmin].mu, particle[lmax-1].mu);
    if (particle[lmin].mu < 0.0)
      fputs("Error: negative magnetic moment encountered.\n", stderr);
    fprintf(stderr, " Jacobian: mu ranges from %le to %le.\n",
	    jacobian[jlmin].mu, jacobian[jlmax-1].mu);

    /* Set up mu bin boundaries */
    nmubins = 25; /* ceil(pow((double)nparts/PPBIN, 1.0/3.0)); */
    if (nmubins < minbins) nmubins = minbins;
    fprintf(stderr, " Using %d bins in mu direction.\n", nmubins);
    mubounds = (double *)malloc((nmubins+1)*sizeof(double));
    mubounds[0] = 0.0;
    if (mu_eqpart) { /* equal particle count bins */
      j = (lmax - lmin)/nmubins;
      for (i=1; i<nmubins; i++)
	mubounds[i] = 0.5*(particle[lmin + i*j - 1].mu + particle[lmin + i*j].mu);
    } else { /* equal width bins */
      delta = (particle[lmax-1].mu - mubounds[0])/nmubins;
      for (i=1; i<nmubins; i++)
	mubounds[i] = mubounds[i-1] + delta;
    }
    mubounds[nmubins] = particle[lmax-1].mu;

    sprintf(outname, "pdist%02d.spl", particle[lmin].sig);
    fp = fopen(outname, "w");
    fprintf(fp, "%d\n", nmubins);

    /* Loop over mu bins */
    for (i=0, binstart=lmin, jbinstart=jlmin; i<nmubins; i++) {
      printf(" Bin %3d: %le < mu < %le (width = %le)\n",
	     i, mubounds[i], mubounds[i+1], mubounds[i+1]-mubounds[i]);
      fprintf(fp, "%.16le\n", mubounds[i+1]);
#ifdef DEBUG
      sprintf(f2name, "mucoarse%02d_%d", particle[lmin].sig, i);
      fp2 = fopen(f2name, "w");
      fprintf(fp2, "%le\t%le\n", mubounds[i], mubounds[i+1]);
#endif

      /* Find last entry in this bin, bounds on momentum, energy */
      pmin = pmax = particle[binstart].pphi;
      emax = particle[binstart].ke;
      for (binstop=binstart; binstop<lmax; binstop++) {
	if (particle[binstop].mu > mubounds[i+1]) break;
	if (particle[binstop].pphi > pmax) pmax = particle[binstop].pphi;
	else if (particle[binstop].pphi < pmin) pmin = particle[binstop].pphi;
	if (particle[binstop].ke > emax) emax = particle[binstop].ke;
      }
      binstop--;
      printf("  1st entry at %d, last at %d (%d total)\n",
	     binstart, binstop, binstop+1-binstart);
      printf("  %le < p_phi < %le\n", pmin, pmax);
      printf("  %le mu < ke < %le\n", cmin, emax);
      fprintf(fp, "%.16le\t%.16le\t%.16le\t%.16le\n", cmin, emax, pmin, pmax);
      for (jbinstop=jbinstart; jbinstop<jlmax; jbinstop++)
	if (jacobian[jbinstop].mu > mubounds[i+1]) break;
      jbinstop--;
      printf("  1st jacobian entry at %d, last at %d (%d total)\n",
	     jbinstart, jbinstop, jbinstop+1-jbinstart);

      /* Compute bin sizes */
      ratio = 1.0 - cmin*mubounds[i]/emax;
      npbins = ceil(sqrt((double)(binstop+1-binstart)/(ratio*PPBIN)));
      nkbins = ceil(ratio*npbins);
      if (npbins < minbins) npbins = minbins;
      if (nkbins < minbins) nkbins = minbins;
      printf("  %d x %d bins...\n", nkbins, npbins);
      pbinwidth = 1.000001*(pmax - pmin)/npbins;
      ebinwidth = 1.000001*(emax - cmin*mubounds[i])/nkbins;

      /* Fill bins */
      fillbins(&binarr, nkbins, cmin, ebinwidth, npbins, pmin, pbinwidth,
	       particle, binstart, binstop);
      /* ko = killorphans(binarr, nkbins, npbins); 
	 printf("  %d orphan(s) removed.\n", ko);
	 printf("  %.2lf ppb\n", (binstop - binstart + 1 - ko)/(double)(npbins*nkbins));
	 kocount += ko; */
      fillbins(&jacarr, nkbins, cmin, ebinwidth, npbins, pmin, pbinwidth,
	       jacobian, jbinstart, jbinstop);
      gaussblur(jacarr, nkbins, npbins, sigma_jac, sigma_jac);

      /* Take ratio */
      for (j=ko=0; j<nkbins; j++) {
	for (k=0; k<npbins; k++) {
	  if (jacarr[j][k] > 0.0)
	    binarr[j][k] /= jacarr[j][k];
	  else {
	    if (binarr[j][k] != 0.0) ko++;
	    binarr[j][k] = 0.0;
	  }
	}
	free(jacarr[j]);
      }
      free(jacarr);
      printf("  %d/%d anomal%s.\n", ko, nkbins*npbins, (ko==1)?"y":"ies");
      gaussblur(binarr, nkbins, npbins, sigma_ptcl, sigma_ptcl);

#ifdef DEBUG
      fprintf(fp2, "%le\t%le\n", pmin, pmax);
      fprintf(fp2, "%le\t%le\n", cmin, emax);
      fprintf(fp2, "%d\t%d\n", nkbins, npbins);
      /* Write bin counts */
      for (j=dtally=0; j<nkbins; j++) {
	for (k=0; k<npbins; k++) {
	  fprintf(fp2, "%le\t", binarr[j][k]);
	  dtally += binarr[j][k];
	}
	fprintf(fp2, "\n");
      }
      fclose(fp2);
      printf("  Total = %lf\n", dtally);
#endif

      /* Allocate, initialize 1D arrays for splining */
      if ((pphi = (double *)malloc((npbins+1) * sizeof(double))) == NULL) {
	fputs("Out of memory\n", stderr);
	exit(1);
      }
      for (k=0; k<npbins; k++)
	pphi[k] = pmin + (k + 0.5)*pbinwidth;
      pphi[npbins] = pmax;
      if ((count = (double *)malloc((npbins+1) * sizeof(double))) == NULL) {
	fputs("Out of memory\n", stderr);
	exit(1);
      }
      count[npbins] = 0.0;
      if ((pspline = (spline **)malloc(nkbins * sizeof(spline *))) == NULL) {
	fputs("Out of memory\n", stderr);
	exit(1);
      }

      /* Construct 1D splines of pphi at each energy */
      fprintf(stderr, "  Creating %d splines...\n", nkbins);
      npcoefs = 5*npbins/8;
      if (npcoefs < 6) npcoefs = 6;
      elvolinv = 1.0/(nparts*(mubounds[i+1]-mubounds[i])*pbinwidth*ebinwidth);
      for (j=0; j<nkbins; j++) { /* For each energy bin */
	for (k=0; k<npbins; k++) /* For each p_phi bin */
	  count[k] = binarr[j][k]*elvolinv;
	free(binarr[j]);

	pspline[j] = createBspline1d(pphi, count, npbins+1, npcoefs, 4, pmin, pmax, 0);
      } /* end loop j */
      free(binarr);  free(count);

      /* Now construct 1D splines of pphi coefficients across energies */
      printf("  Splining %d splines...\n", nkbins);
      if ((count = (double *)malloc((nkbins+1) * sizeof(double))) == NULL) {
	fputs("Out of memory\n", stderr);
	exit(1);
      }
      if ((ke = (double *)malloc(nkbins * sizeof(double))) == NULL) {
	fputs("Out of memory\n", stderr);
	exit(1);
      }
      for (k=0; k<nkbins; k++)
	ke[k] = cmin*mubounds[i] + (k + 0.5)*ebinwidth;
      nkcoefs = 5*nkbins/8;
      if (nkcoefs < 6) nkcoefs = 6;
      fprintf(fp, "%d\t%d\n", npcoefs, nkcoefs);
      for (j=0; j<npcoefs; j++) { /* For each coefficient */
	for (k=0; k<nkbins; k++) /* For each ke bin */
	  count[k] = pspline[k]->coefs[j];

	kspline = createBspline1d(ke, count, nkbins, nkcoefs, 4, cmin*mubounds[i], emax, 0);

	for (k=0; k<nkcoefs; k++)
	  fprintf(fp, "%.16le\t", kspline->coefs[k]);
	fprintf(fp, "\n");

	splinedelete(kspline); free(kspline);
      } /* end loop j */

      /* Free up bin storage */
      free(pphi);  free(count);  free(ke);
      for (j=0; j<nkbins; j++) { splinedelete(pspline[j]); free(pspline[j]);}
      free(pspline);

      binstart = binstop + 1;  jbinstart = jbinstop + 1;
    } /* end loop i */

    fclose(fp);
    free(mubounds);

    lmin = lmax;  lmax = nparts;
    jlmin = jlmax; jlmax = njac;
  } while (lmin < lmax);

  printf("%d / %d particles (%.6f %%) deleted.\n", kocount, nparts,
	  (100.0*kocount)/nparts);

  free(particle); free(jacobian);

  return 0;
}

/******************************************************************************/
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
  fprintf(stderr, "%ld species, %ld particles\n", (long)len1, (long)len0);
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

/******************************************************************************/
int mcompar(const void *f1, const void *f2)
{
  double mu1, mu2;
  int    s1 = ((vvect *)f1)->sig, s2 = ((vvect *)f2)->sig;

  if (s1==s2) {
    mu1 = ((vvect *)f1)->mu;  mu2 = ((vvect *)f2)->mu;
    if (mu1 < mu2) return -1;
    if (mu1 > mu2) return 1;
    return 0;
  }

  return s1-s2;
}

/******************************************************************************/
void fillbins(double ***arr, int nk, double cmin, double ebinwidth,
	      int np, double pmin, double pbinwidth,
	      vvect *part, int ji, int jf)
{
  int j, ik, ip;

  /* Allocate storage, initialize */
  if ((*arr = (double **)malloc(nk * sizeof(double *))) == NULL) {
    fputs("Out of memory in fillbins\n", stderr);
    exit(1);
  }
  for (j=0; j<nk; j++) {
    if (((*arr)[j] = (double *)malloc(np * sizeof(double))) == NULL) {
      fputs("Out of memory in fillbins\n", stderr);
      exit(1);
    }
    memset((*arr)[j], 0, np * sizeof(double));
  }

  /* Loop through particles, assigning to bins */
  for (j=ji; j<=jf; j++) {
    ik = (int)((part[j].ke - cmin*part[j].mu)/ebinwidth);
    if ((ik>=0) && (ik<nk)) {
      ip = (int)((part[j].pphi - pmin)/pbinwidth);
      if ((ip>=0) && (ip<np)) (*arr)[ik][ip] += part[j].wt;
    }
  }
}

/******************************************************************************/
/* Filter some noise by removing isolated particles, presumed lost.           */
int killorphans(int **barr, const int m, const int n)
{
  int row, col, tally=0;

  fputs("  Removing orphans...\n", stderr);

  /* Array interior */
  for (row=1; row<m-1; row++)
    for (col=1; col<n-1; col++)
      if (barr[row][col] == 1)
	if ((barr[row-1][col] == 0) && (barr[row][col-1] == 0) &&
	    (barr[row][col+1] == 0) && (barr[row+1][col] == 0)) {
	  barr[row][col] = 0; ++tally;
	}

  /* First & last rows */
  for (col=1; col<n-1; col++) {
    if (barr[0][col] == 1)
      if ((barr[0][col-1] == 0) && (barr[0][col+1] == 0) &&
	  (barr[1][col] == 0)) {
	barr[0][col] = 0; ++tally;
      }
    if (barr[m-1][col] == 1)
      if ((barr[m-1][col-1] == 0) && (barr[m-1][col+1] == 0) &&
	  (barr[m-2][col] == 0)) {
	barr[m-1][col] = 0; ++tally;
      }
  }

  /* First & last cols */
  for (row=1; row<m-1; row++) {
    if (barr[row][0] == 1)
      if ((barr[row-1][0] == 0) &&
	  (barr[row][1] == 0) && (barr[row+1][0] == 0)) {
	barr[row][0] = 0; ++tally;
      }
    if (barr[row][n-1] == 1)
      if ((barr[row-1][n-1] == 0) && (barr[row][n-2] == 0) &&
	  (barr[row+1][n-1] == 0)) {
	barr[row][n-1] = 0; ++tally;
      }
  }

  /* Corners */
  if ((barr[0][0] == 1) && (barr[0][1] == 0) && (barr[1][0] == 0)) {
    barr[0][0] = 0; ++tally;
  }
  if ((barr[0][n-1] == 1) && (barr[0][n-2] == 0) && (barr[1][n-1] == 0)) {
    barr[0][n-1] = 0; ++tally;
  }
  if ((barr[m-1][0] == 1) && (barr[m-2][0] == 0) && (barr[m-1][1] == 0)) {
    barr[m-1][0] = 0; ++tally;
  }
  if ((barr[m-1][n-1] == 1) && (barr[m-2][n-1] == 0) && (barr[m-1][n-2] == 0)) {
    barr[m-1][n-1] = 0; ++tally;
  }

  return tally;
}

/******************************************************************************/
int threshold(int **barr, const int nrows, const int ncols, const int mincount)
{
  int row, col, tally=0;

  for (row=0; row<nrows; row++)
    for (col=0; col<ncols; col++)
      if (barr[row][col] < mincount) {
	tally += barr[row][col];
	barr[row][col] = 0;
      }

  return tally;
}

/******************************************************************************/
void gaussblur(double **arr, const int nrows, const int ncols,
		const double sigma_x, const double sigma_y)
{
  double **work;
  double  *kernel;
  double   sum, twossqinv;
  int      kernel_size, ik, x, irow, icol, offset, cstart, cstop, kstart, kstop;

  fprintf(stderr, "  gaussblur called; sigma = %lf, %lf\n", sigma_x, sigma_y);

  /* Set up temporary workspace */
  if ((work = (double **)malloc(nrows * sizeof(double *))) == NULL) {
    fputs("Error: out of memory in gaussblur.\n", stderr);
    return;
  }
  for (irow=0; irow<nrows; irow++)
    if ((work[irow] = (double *)malloc(ncols * sizeof(double))) == NULL) {
      fputs("Error: out of memory in gaussblur.\n", stderr);
      return;
    }

  /* Calculate x kernel size, allocate storage */
  kernel_size = 1 + 2*(int)(3.0*sigma_x); /* ensures odd number, so center gets highest weight */
  if (kernel_size < 2) {
    for (irow=0; irow<nrows; irow++)
      memcpy(work[irow], arr[irow], ncols*sizeof(double));
  } else {
    printf("  Blurring rows using kernel size %d...\n", kernel_size);
    if ((kernel = (double *)malloc(kernel_size * sizeof(double))) == NULL) {
      fputs("Error: out of memory in gaussblur.\n", stderr);
      return;
    }

    /* Calculate, normalize x kernel */
    twossqinv = 0.5/(sigma_x*sigma_x);
    /* norm = sqrt(twossqinv / 3.1415926535897932384626); */
    for (ik=0, sum=0.0; ik<kernel_size; ik++) {
      x = ik - kernel_size/2;
      kernel[ik] = exp(-twossqinv*x*x);
      sum += kernel[ik];
    }
    for (ik=0; ik<kernel_size; ik++) kernel[ik] /= sum;

    /* Apply convolution to each row */
    for (irow=0; irow<nrows; irow++)
      for (icol=0; icol<ncols; icol++) {

	offset = cstart = icol - kernel_size/2;
	if (cstart < 0) cstart = 0;
	kstart = cstart - offset;
	cstop = icol + kernel_size/2;
	if (cstop > ncols-1) cstop = ncols - 1;
	kstop = kstart + (cstop - cstart);

	work[irow][icol] = 0.0;
	for (ik=0; ik<kstart; ik++)
	  work[irow][icol] += kernel[ik] * arr[irow][0];
	for (; ik<=kstop; ik++)
	  work[irow][icol] += kernel[ik] * arr[irow][ik + offset];
	for (; ik<kernel_size; ik++)
	  work[irow][icol] += kernel[ik] * arr[irow][ncols - 1];

      } /* end loop icol */

    free(kernel);
  } /* end else (sigma_x is positive) */

  /* Calculate y kernel size, allocate storage */
  kernel_size = 1 + 2*(int)(3.0*sigma_y); /* ensures odd number, so center gets highest weight */
  if (kernel_size < 2) {
    for (irow=0; irow<nrows; irow++)
      memcpy(arr[irow], work[irow], ncols*sizeof(double));
  } else {
    printf("  Blurring columns using kernel size %d...\n", kernel_size);
    if ((kernel = (double *)malloc(kernel_size * sizeof(double))) == NULL) {
      fputs("Error: out of memory in gaussblur.\n", stderr);
      return;
    }

    /* Calculate, normalize y kernel */
    twossqinv = 0.5/(sigma_y*sigma_y);
    for (ik=0, sum=0.0; ik<kernel_size; ik++) {
      x = ik - kernel_size/2;
      kernel[ik] = exp(-twossqinv*x*x);
      sum += kernel[ik];
    }
    for (ik=0; ik<kernel_size; ik++) kernel[ik] /= sum;

    /* Now apply convolution to each column */
    for (icol=0; icol<ncols; icol++)
      for (irow=0; irow<nrows; irow++) {

	offset = cstart = irow - kernel_size/2;
	if (cstart < 0) cstart = 0;
	kstart = cstart - offset;
	cstop = irow + kernel_size/2;
	if (cstop > nrows-1) cstop = nrows - 1;
	kstop = kstart + (cstop - cstart);

	arr[irow][icol] = 0.0;
	for (ik=0; ik<kstart; ik++)
	  arr[irow][icol] += kernel[ik] * work[0][icol];
	for (; ik<=kstop; ik++)
	  arr[irow][icol] += kernel[ik] * work[ik + offset][icol];
	for (; ik<kernel_size; ik++)
	  arr[irow][icol] += kernel[ik] * work[nrows - 1][icol];
      } /* end loop irow */

    free(kernel);
  } /* end else (sigma_y is positive) */

  for (irow=0; irow<nrows; irow++) free(work[irow]);
  free(work);
  fputs("  Blur applied.\n", stderr);
}

/******************************************************************************/
int getlbound(vvect *pvec, const int n)
{
  int j;

  for (j=0; j<n; j++)
    if (pvec[j].sig > 0) break;

  return j;
}
