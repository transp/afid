#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "spline_interface.h"
#include "particleIO.h"

/* Function prototypes */
void fillbins(double ***arr, int nk, double cmin, double ew,
	      int np, double pmin, double pw,
	      vvect *part, int ji, int jf);
int killorphans(double **, const int, const int);
int threshold(int **, const int, const int, const int);
void gaussblur(double **, const int, const int, const double, const double, const int);
int getlbound(vvect *, const int);

/******************************************************************************/
int main(int argc, char *argv[])
{
  FILE   *fp, *fp2;
  char outname[32], inname[96], f2name[64];
  vvect  *particle, *jacobian;
  spline **pspline;
  spline *kspline;
  double *ke, *pphi, *count;
  double *mubounds;
  double **binarr, **jacarr;
  double  sigma_ptcl = 1.334, sigma_jac = 1.75;
  double  delta, ratio, pmin, pmax, cmin, emax;
  double  pbinwidth, ebinwidth, elvolinv, dtally, minlog, cbrat=0.625, rsq, rsqtol=0.985;
  int     nparts, njac, nmubins=25, i, j, k, nkbins, npbins, npcoefs, nkcoefs;
  int     binstart, binstop, jbinstart, jbinstop, ko, kocount=0, lbound, jlbound;
  int     lmin=0, lmax, jlmin=0, jlmax, mu_eqpart=1, ppbin=40, minbins=7, debug_flag=0;
  int     logbins=0, nzero, koflag=0, bincount=0, mincoefs=6, csorted, jsorted;

  /* If there are no arguments, print a usage message. */
  if (argc==1) {
    fputs("Usage: fitjac <fileroot> [options]\n", stderr);
    fputs("\n\tBins particle distribution data from <fileroot>_condensed.cdf,\n", stderr);
    fputs("\tnormalizes according to data in <fileroot>_jacobian.cdf, and fits\n", stderr);
    fputs("\tsplines in E and P_phi to a series of mu cross-sections.\n", stderr);
    fputs("\nOptions:\n", stderr);
    fprintf(stderr, "\t-log: fit splines to ln(f) rather than f itself.\n");
    fprintf(stderr, "\t-minbins <min. # of bins in any direction> (default %d)\n", minbins);
    fprintf(stderr, "\t-ppbin <particles/bin target> (default %d)\n", ppbin);
    fprintf(stderr, "\t-nmu <number of mu bins> (default %d), 0 to calculate based on ppbin\n", nmubins);
    fprintf(stderr, "\t-mueqw: specify equal-width mu bins (default: equal particle count)\n");
    fprintf(stderr, "\t-ii: ignore isolated single-particle bins\n");
    fprintf(stderr, "\t-sigma_p <particle bin gauss blur radius> (default %.3lf)\n", sigma_ptcl);
    fprintf(stderr, "\t-sigma_j <Jacobian bin gauss blur radius> (default %.3lf)\n", sigma_jac);
    fprintf(stderr, "\t-crat <ratio of coef count to bin count> (default %.3lf)\n", cbrat);
    fprintf(stderr, "\t-mincoefs <min. # of spline coefs in any direction> (default %d)\n", mincoefs);
    fprintf(stderr, "\t-debug: generate debugging information, including mucoarse files.\n");
    fputs("\n", stderr);
    exit(0);
  }

  /* Parse command-line options */
  for (i=2; i<argc; i++) {
    if (!strncmp(argv[i], "-log", 4)) {
      logbins = 1;
      fputs("Fitting log of distribution.\n", stderr);
      continue;
    }
    if (!strncmp(argv[i], "-minbins", 8)) {
      minbins = atoi(argv[++i]);
      if (minbins < 7) minbins = 7;
      fprintf(stderr, "Setting min. bin count to %d.\n", minbins);
      continue;
    }
    if (!strncmp(argv[i], "-ppbin", 6)) {
      ppbin = atoi(argv[++i]);
      if (ppbin < 2) ppbin = 2;
      fprintf(stderr, "Targeting %d particles/bin.\n", ppbin);
      continue;
    }
    if (!strncmp(argv[i], "-nmu", 4)) {
      nmubins = atoi(argv[++i]);
      continue;
    }
    if (!strncmp(argv[i], "-mueqw", 6)) {
      mu_eqpart = 0;
      fputs("Using equally spaced mu bins.\n", stderr);
      continue;
    }
    if (!strncmp(argv[i], "-ii", 3)) {
      koflag = 1;
      fputs("Ignoring isolated single-occupancy bins.\n", stderr);
      continue;
    }
    if (!strncmp(argv[i], "-sigma_p", 8)) {
      sigma_ptcl = atof(argv[++i]);
      if (sigma_ptcl < 0.0) sigma_ptcl = 0.0;
      fprintf(stderr, "Setting particle blurring radius to %.6lf.\n", sigma_ptcl);
      continue;
    }
    if (!strncmp(argv[i], "-sigma_j", 8)) {
      sigma_jac = atof(argv[++i]);
      if (sigma_jac < 0.0) sigma_jac = 0.0;
      fprintf(stderr, "Setting Jacobian blurring radius to %.6lf.\n", sigma_jac);
      continue;
    }
    if (!strncmp(argv[i], "-crat", 5)) {
      cbrat = atof(argv[++i]);
      fprintf(stderr, "Setting coef-to-bin ratio to %.6lf.\n", cbrat);
      continue;
    }
    if (!strncmp(argv[i], "-mincoefs", 9)) {
      mincoefs = atoi(argv[++i]);
      fprintf(stderr, "Setting min. coef count to %d.\n", mincoefs);
      continue;
    }
    if (!strncmp(argv[i], "-debug", 6)) {
      debug_flag = 1;
      fputs("Generating debugging info.\n", stderr);
      continue;
    }
  }

  /* Read all particle data from file */
  sprintf(inname, "%s_condensed.cdf", argv[1]);
  readParticleData(inname, &particle, &nparts, &csorted, 0);
  printf("%d %ssorted particles read; last = %le  %le  %le  %d  %.3le\n", nparts,
	 (csorted) ? "" : "un",
	 particle[nparts-1].pphi, particle[nparts-1].mu, particle[nparts-1].ke,
	 particle[nparts-1].sig, particle[nparts-1].wt);

  /* Read all jacobian data from file */
  sprintf(inname, "%s_jacobian.cdf", argv[1]);
  readParticleData(inname, &jacobian, &njac, &jsorted, 0);
  printf("%d %ssorted particles read; last = %le  %le  %le  %d  %.3le\n", njac,
	 (jsorted) ? "" : "un",
	 jacobian[njac-1].pphi, jacobian[njac-1].mu, jacobian[njac-1].ke,
	 jacobian[njac-1].sig, jacobian[njac-1].wt);

  /* Calculate cmin (minimum ratio of ke to mu), emax */
  cmin = particle[0].ke / particle[0].mu;  emax = particle[0].ke;
  for (i=1; i<nparts; i++) {
    if (particle[i].ke > emax) emax = particle[i].ke;
    ratio = particle[i].ke / particle[i].mu;
    if (ratio < cmin) cmin = ratio;
  }
  printf("cmin = %lf; emax = %le\n", cmin, emax);

  /* Sort the particles by sign of v||, increasing magnetic moment mu */
  if (!csorted) {
    fputs("Sorting by lambda, magnetic moment...\n", stderr);
    qsort(particle, nparts, sizeof(vvect), mcompar);
  }
  lmax = lbound = getlbound(particle, nparts);
  printf("%d condensed particles have l<0 (%.1lf %%).\n", lbound, 100.0*lbound/nparts);
  if (!jsorted) {
    fputs("Sorting Jacobian...\n", stderr);
    qsort(jacobian, njac, sizeof(vvect), mcompar);
  }
  jlmax = jlbound = getlbound(jacobian, njac);
  printf("%d Jacobian particles have l<0 (%.1lf %%).\n", jlbound, 100.0*jlbound/njac);

  do { /* Once for lambda<0, once for lambda>0 */
    printf("sgn=%d: mu ranges from %le to %le.\n", particle[lmin].sig,
	   particle[lmin].mu, particle[lmax-1].mu);
    if (particle[lmin].mu < 0.0)
      fputs("Error: negative magnetic moment encountered.\n", stderr);
    printf("Jacobian: mu ranges from %le to %le.\n",
	   jacobian[jlmin].mu, jacobian[jlmax-1].mu);

    /* Set up mu bin boundaries */
    if (nmubins < 1) nmubins = ceil(pow((double)nparts/ppbin, 1.0/3.0));
    if (nmubins < minbins) nmubins = minbins;
    printf("Using %d bins in mu direction.\n", nmubins);
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
    if (logbins) fputs("Log\n", fp);
    fprintf(fp, "%d\n", nmubins);

    /* Loop over mu bins */
    for (i=0, binstart=lmin, jbinstart=jlmin; i<nmubins; i++) {
      printf(" Bin %3d: %le < mu < %le (width = %le)\n",
	     i, mubounds[i], mubounds[i+1], mubounds[i+1]-mubounds[i]);
      fprintf(fp, "%.16le\n", mubounds[i+1]);

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
      npbins = ceil(sqrt((double)(binstop+1-binstart)/(ratio*ppbin)));
      nkbins = ceil(ratio*npbins);
      if (npbins < minbins) npbins = minbins;
      if (nkbins < minbins) nkbins = minbins;
      printf("  %d x %d bins...\n", nkbins, npbins);
      pbinwidth = 1.000001*(pmax - pmin)/npbins;
      ebinwidth = 1.000001*(emax - cmin*mubounds[i])/nkbins;

      /* Fill bins */
      fillbins(&binarr, nkbins, cmin, ebinwidth, npbins, pmin, pbinwidth,
	       particle, binstart, binstop);
      bincount += nkbins*npbins;
      if (koflag) {
	kocount += (ko = killorphans(binarr, nkbins, npbins));
	if (ko || debug_flag)
	  printf("  %d isolated bin%s zeroed.\n", ko, (ko>1)?"s":"");
      } else ko=0;
      printf("  Avg %.2lf particles/bin\n", (binstop - binstart + 1 - ko)/(double)(npbins*nkbins));
      fillbins(&jacarr, nkbins, cmin, ebinwidth, npbins, pmin, pbinwidth,
	       jacobian, jbinstart, jbinstop);
      gaussblur(jacarr, nkbins, npbins, sigma_jac, sigma_jac, debug_flag);

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
      if (ko) printf("  %d/%d anomal%s.\n", ko, nkbins*npbins, (ko==1)?"y":"ies");
      gaussblur(binarr, nkbins, npbins, sigma_ptcl, sigma_ptcl, debug_flag);

      if (debug_flag) {
	sprintf(f2name, "mucoarse%02d_%d", particle[lmin].sig, i);
	if ((fp2 = fopen(f2name, "w")) == NULL)
	  fprintf(stderr, "Error: could not create file %s for writing.\n", f2name);
	else {
	  /* Write header info to mucoarse file */
	  fprintf(fp2, "%le\t%le\n", mubounds[i], mubounds[i+1]);
	  fprintf(fp2, "%le\t%le\n", pmin, pmax);
	  fprintf(fp2, "%le\t%le\n", cmin, emax);
	  fprintf(fp2, "%d\t%d\n", nkbins, npbins);

	  /* Write bin counts */
	  dtally = 0.0;
	  for (j=0; j<nkbins; j++) {
	    for (k=0; k<npbins; k++) {
	      fprintf(fp2, "%le\t", binarr[j][k]);
	      dtally += binarr[j][k];
	    }
	    fprintf(fp2, "\n");
	  }
	  fclose(fp2);
	  printf("  Total = %le\n", dtally);
	}
      } // end if debug_flag

      // Renormalize
      elvolinv = 1.0/(nparts*(mubounds[i+1]-mubounds[i])*pbinwidth*ebinwidth);
      for (j=0; j<nkbins; j++)   /* For each energy bin */
	for (k=0; k<npbins; k++) /* For each p_phi bin */
	  binarr[j][k] *= elvolinv;

      if (logbins) { // Compute natural log
	minlog = 1.0e+299;
	for (nzero=j=0; j<nkbins; j++)
	  for (k=0; k<npbins; k++)
	    if (binarr[j][k] > 0.0) {
	      binarr[j][k] = log(binarr[j][k]);
	      if (binarr[j][k] < minlog) minlog = binarr[j][k];
	    } else ++nzero; // end loop k, j
	if (debug_flag) {
	  printf("  minlog = %lf\n", minlog);
	  printf("  %d zero val(s) to be replaced.\n", nzero);
	}
	if (nzero) {
	  minlog -= 6.907755;
	  for (j=0; j<nkbins; j++)
	    for (k=0; k<npbins; k++)
	      if (binarr[j][k] <= 0.0) binarr[j][k] = minlog;
	}
      } // end if logbins

      /* Allocate, initialize 1D arrays for splining */
      if ((pphi = (double *)malloc((npbins+1) * sizeof(double))) == NULL) {
	fputs("Out of memory\n", stderr);
	exit(1);
      }
      for (k=0; k<npbins; k++)
	pphi[k] = pmin + (k + 0.5)*pbinwidth;
      pphi[npbins] = (logbins) ? pmin : pmax;
      if ((count = (double *)malloc((npbins+1) * sizeof(double))) == NULL) {
	fputs("Out of memory\n", stderr);
	exit(1);
      }
      if ((pspline = (spline **)malloc(nkbins * sizeof(spline *))) == NULL) {
	fputs("Out of memory\n", stderr);
	exit(1);
      }

      /* Construct 1D splines of pphi at each energy */
      npcoefs = ceil(cbrat*npbins);
      if (npcoefs < mincoefs) npcoefs = mincoefs;
      fprintf(stderr, "  mu bin %3d_%02d: Creating %d splines of %d coefs...\n",
	      i, particle[lmin].sig, nkbins, npcoefs);
      for (j=0; j<nkbins; j++) { /* For each energy bin */
	for (k=0; k<npbins; k++) /* For each p_phi bin */
	  count[k] = binarr[j][k];
	count[npbins] = (logbins) ? count[0] : 0.0;
	free(binarr[j]);

	pspline[j] = createBspline1d(pphi, count, npbins +1, npcoefs, 4,
				     pmin, pmax, &rsq, 0);
	if (debug_flag) printf("  R^2 = %lf\n", rsq);
	if (rsq < rsqtol) fprintf(stderr, "Warning: poor fit quality: R^2 = %lf\n", rsq);
      } /* end loop j */
      free(binarr);  free(count);

      /* Now construct 1D splines of pphi coefficients across energies */
      if ((count = (double *)malloc(nkbins * sizeof(double))) == NULL) {
	fputs("Out of memory\n", stderr);
	exit(1);
      }
      if ((ke = (double *)malloc(nkbins * sizeof(double))) == NULL) {
	fputs("Out of memory\n", stderr);
	exit(1);
      }
      for (k=0; k<nkbins; k++)
	ke[k] = cmin*mubounds[i] + (k + 0.5)*ebinwidth;
      nkcoefs = ceil(cbrat*nkbins);
      if (nkcoefs < mincoefs) nkcoefs = mincoefs;
      if (debug_flag)
	fprintf(stderr, "  Splining %d splines with %d coefs...\n", nkbins, nkcoefs);
      fprintf(fp, "%d\t%d\n", npcoefs, nkcoefs);
      for (j=0; j<npcoefs; j++) { /* For each coefficient */
	for (k=0; k<nkbins; k++) /* For each ke bin */
	  count[k] = pspline[k]->coefs[j];

	kspline = createBspline1d(ke, count, nkbins, nkcoefs, 4,
				  cmin*mubounds[i], emax, &rsq, 0);
	printf("  R^2[%d] = %lf\n", j, rsq);
	if (rsq < rsqtol) fputs("Warning: poor fit quality!\n", stderr);

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

  if (koflag) printf("%d / %d bins (%.6f %%) zeroed.\n", kocount, bincount,
		     (100.0*kocount)/bincount);

  free(particle); free(jacobian);
  fputs("Done.\n", stderr);
  return 0;
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
int killorphans(double **barr, const int m, const int n)
{
  int row, col, tally=0;

  /* Array interior */
  for (row=1; row<m-1; row++)
    for (col=1; col<n-1; col++)
      if (barr[row][col] > 0.0)
	if ((barr[row-1][col] == 0.0) && (barr[row][col-1] == 0.0) &&
	    (barr[row][col+1] == 0.0) && (barr[row+1][col] == 0.0)) {
	  barr[row][col] = 0.0; ++tally;
	}

  /* First & last rows */
  for (col=1; col<n-1; col++) {
    if (barr[0][col] > 0.0)
      if ((barr[0][col-1] == 0.0) && (barr[0][col+1] == 0.0) &&
	  (barr[1][col] == 0.0)) {
	barr[0][col] = 0.0; ++tally;
      }
    if (barr[m-1][col] > 0.0)
      if ((barr[m-1][col-1] == 0.0) && (barr[m-1][col+1] == 0.0) &&
	  (barr[m-2][col] == 0.0)) {
	barr[m-1][col] = 0.0; ++tally;
      }
  }

  /* First & last cols */
  for (row=1; row<m-1; row++) {
    if (barr[row][0] > 0.0)
      if ((barr[row-1][0] == 0.0) &&
	  (barr[row][1] == 0.0) && (barr[row+1][0] == 0.0)) {
	barr[row][0] = 0.0; ++tally;
      }
    if (barr[row][n-1] > 0.0)
      if ((barr[row-1][n-1] == 0.0) && (barr[row][n-2] == 0.0) &&
	  (barr[row+1][n-1] == 0.0)) {
	barr[row][n-1] = 0.0; ++tally;
      }
  }

  /* Corners */
  if ((barr[0][0] > 0.0) && (barr[0][1] == 0.0) && (barr[1][0] == 0.0)) {
    barr[0][0] = 0.0; ++tally;
  }
  if ((barr[0][n-1] > 0.0) && (barr[0][n-2] == 0.0) && (barr[1][n-1] == 0.0)) {
    barr[0][n-1] = 0.0; ++tally;
  }
  if ((barr[m-1][0] > 0.0) && (barr[m-2][0] == 0.0) && (barr[m-1][1] == 0.0)) {
    barr[m-1][0] = 0.0; ++tally;
  }
  if ((barr[m-1][n-1] > 0.0) && (barr[m-2][n-1] == 0.0) && (barr[m-1][n-2] == 0.0)) {
    barr[m-1][n-1] = 0.0; ++tally;
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
	       const double sigma_x, const double sigma_y, const int df)
{
  double **work;
  double  *kernel;
  double   sum, twossqinv;
  int      kernel_size, ik, x, irow, icol, offset, cstart, cstop, kstart, kstop;

  if (df) fprintf(stderr, "  gaussblur called; sigma = %lf, %lf\n", sigma_x, sigma_y);

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
    if (df) printf("  Blurring rows using kernel size %d...\n", kernel_size);
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
    if (df) printf("  Blurring columns using kernel size %d...\n", kernel_size);
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
  if (df) fputs("  Blur applied.\n", stderr);
}

/******************************************************************************/
int getlbound(vvect *pvec, const int n)
{
  int j;

  for (j=0; j<n; j++)
    if (pvec[j].sig > 0) break;

  return j;
}
