#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "particleSplines.h"
#ifdef PARTICLE_MPI
#include "mpi.h"
#endif

#define PDF_NEAREST_NEIGHBOR
#ifdef PDF_NEAREST_NEIGHBOR
#define PDFFN pdfnn
#define PDFDR pdfnnder
#else
#define PDFFN pdflin
#define PDFDR pdflinder
#endif

spline3d *psplinedata, *nsplinedata;

/* Fortran interface routines */

/* Read particle spline info from file to intialize data structure; */
/*  return value is number of 2D splines read.                      */
int pspline_init_(int *n2d, int *p2d)
{
  double posint, negint;
  const int nsplint=192, nmuint=64;
  const char pname[] = "pdist01.spl", nname[] = "pdist-1.spl";

#ifdef PDF_NEAREST_NEIGHBOR
  fputs("Setting up splines for nearest-neighbor mu interpolation.\n", stderr);
#else
  fputs("Setting up splines for linear mu interpolation.\n", stderr);
#endif

  /* Read data from files */
  nsplinedata = readspline(nname);  psplinedata = readspline(pname);

  /* Compute appropriate norm so that f integrates to unity */
  posint = integratedV(psplinedata, psplinedata->pmin, psplinedata->pmax, nsplint,
	  psplinedata->mubounds[0], psplinedata->mubounds[psplinedata->nmubins], nmuint,
	  psplinedata->cmin*psplinedata->mubounds[0], psplinedata->emax, nsplint);
  negint = integratedV(nsplinedata, nsplinedata->pmin, nsplinedata->pmax, nsplint,
	  nsplinedata->mubounds[0], nsplinedata->mubounds[nsplinedata->nmubins], nmuint,
	  nsplinedata->cmin*nsplinedata->mubounds[0], nsplinedata->emax, nsplint);
  psplinedata->norm = nsplinedata->norm = 1.0/(posint + negint);

  *n2d = nsplinedata->nmubins;  *p2d = psplinedata->nmubins;

  return 0;
}

/******************************************************************************/
/* Return bounds of constants-of-motion-space cube containing the particle */
/*  distribution function. */
int getpspline3bounds_(double *pmin, double *pmax, double *mumin, double *mumax,
		       double *cmin, double *emax)
{
  *pmin = (psplinedata->pmin < nsplinedata->pmin) ?
    psplinedata->pmin : nsplinedata->pmin;
  *pmax = (psplinedata->pmax > nsplinedata->pmax) ?
    psplinedata->pmax : nsplinedata->pmax;
  *mumin = (psplinedata->mubounds[0] < nsplinedata->mubounds[0]) ?
    psplinedata->mubounds[0] : nsplinedata->mubounds[0];
  *mumax = (psplinedata->mubounds[psplinedata->nmubins] >
	    nsplinedata->mubounds[nsplinedata->nmubins]) ?
    psplinedata->mubounds[psplinedata->nmubins] :
    nsplinedata->mubounds[nsplinedata->nmubins];
  *cmin = (psplinedata->cmin < nsplinedata->cmin) ?
    psplinedata->cmin : nsplinedata->cmin;
  *emax = (psplinedata->emax > nsplinedata->emax) ?
    psplinedata->emax : nsplinedata->emax;

  return 0;
}

/******************************************************************************/
/* Return bounds of constants-of-motion-space rectangle containing particle */
/*  distribution function at specified mu, sign(v) */
int getpspline2bounds_(double *mu, int *sgnv,
		       double *pmin, double *pmax, double *cmin, double *emax)
{
  spline3d *data;
  int       mbin;

  data = (*sgnv > 0) ? psplinedata : nsplinedata;

  /* Find the mu bin containing this point */
  for (mbin=0; mbin<data->nmubins; mbin++)
    if (*mu < data->mubounds[mbin+1]) break;
  if (mbin == data->nmubins) {
    fputs("mu out of range in getpspline2bounds\n", stderr);
    return 0;
  }

  /* Look up values */
  *pmin = data->pespline[mbin].pmin;
  *pmax = data->pespline[mbin].pmax;
  *cmin = data->pespline[mbin].cmin;
  *emax = data->pespline[mbin].emax;

  return 0;
}

/******************************************************************************/
/* Deallocate memory for particle spline data */
int pspline_free_()
{
  freespline3d(psplinedata);
  freespline3d(nsplinedata);
  return 0;
}

/******************************************************************************/
/* Evaluate f(p_phi, mu, E) by 2D cubic B-spline interpolation in P,E */
int getpdf_(double *pphi, double *mu, double *ke, int *sgnv, double *f)
{
  if (*sgnv > 0)
    *f = PDFFN(*pphi, *mu, *ke, psplinedata);
  else
    *f = PDFFN(*pphi, *mu, *ke, nsplinedata);

  return 0;
}

/******************************************************************************/
/* Evaluate f(p_phi, mu, E) and its 1st derivatives with respect to P and E */
/*  by 2D cubic B-spline interpolation in P,E */
int getpdfd_(double *pphi, double *mu, double *ke, int *sgnv, double *f,
	    double *dfdp, double *dfdE)
{
  double drv[2];

  if (*sgnv > 0)
    *f = PDFDR(*pphi, *mu, *ke, drv, psplinedata);
  else
    *f = PDFDR(*pphi, *mu, *ke, drv, nsplinedata);

  *dfdp = drv[0];  *dfdE = drv[1];

  return 0;
}

#ifndef PARTICLE_MPI
/******************************************************************************/
/* Particle spline utilities */
/* Serial version of spline coefficient reader */
spline3d *readspline(const char *fname)
{
  FILE     *fp;
  spline3d *data;
  char     buf[64];
  double   emin;
  int      mbin, pcoef, kcoef, npcoefs, nkcoefs;

  /* Allocate storage */
  data = (spline3d *)malloc(sizeof(spline3d));
  data->norm = data->pmin = 1.0;  data->pmax = data->emax = -1.0;  data->cmin=1.0e+9;
  data->logbins = 0;

  /* Open input file */
  if ((fp = fopen(fname, "r")) == NULL) {
    fprintf(stderr, "Error: could not open file %s for reading.\n", fname);
    exit(1);
  }

  /* Check for log flag */
  fgets(buf, 64, fp);
  if ((*buf == 'L') || (*buf == 'l')) {
    data->logbins = 1;
    fprintf(stderr, "%s: bins are logarithmic.\n", fname);
  } else rewind(fp);

  /* Read data */
  fscanf(fp, "%d", &data->nmubins);
  if (data->nmubins < 1) exit(1);
  if ((data->mubounds = (double *)malloc((data->nmubins+1) * sizeof(double))) == NULL) {
    fputs("Out of memory in readspline.\n", stderr);
    exit(2);
  }
  data->mubounds[0] = 0.0;
  if ((data->pespline = (spline2d *)malloc(data->nmubins * sizeof(spline2d))) == NULL) {
    fputs("Out of memory in readspline.\n", stderr);
    exit(2);
  }
  for (mbin=0; mbin<data->nmubins; mbin++) {
    fscanf(fp, "%le", &data->mubounds[mbin+1]);
    if (data->mubounds[mbin] >= data->mubounds[mbin+1]) exit(1);
    fscanf(fp, "%le\t%le\t%le\t%le",
	   &data->pespline[mbin].cmin, &data->pespline[mbin].emax,
	   &data->pespline[mbin].pmin, &data->pespline[mbin].pmax);
    if (data->pespline[mbin].pmin < data->pmin) data->pmin = data->pespline[mbin].pmin;
    if (data->pespline[mbin].pmax > data->pmax) data->pmax = data->pespline[mbin].pmax;
    if (data->pespline[mbin].cmin < data->cmin) data->cmin = data->pespline[mbin].cmin;
    if (data->pespline[mbin].emax > data->emax) data->emax = data->pespline[mbin].emax;
    if (data->pespline[mbin].pmin >= data->pespline[mbin].pmax) exit(1);
    emin = data->pespline[mbin].cmin * data->mubounds[mbin];
    fscanf(fp, "%d\t%d", &npcoefs, &nkcoefs);
    if ((npcoefs < 4) || (nkcoefs < 4)) exit(1);
    data->pespline[mbin].npcoefs = npcoefs;
    data->pespline[mbin].nkcoefs = nkcoefs;
    splinealloc(&data->pespline[mbin].kspline, npcoefs, nkcoefs, 4,
		emin, data->pespline[mbin].emax);
    for (pcoef=0; pcoef<npcoefs; pcoef++)
      for (kcoef=0; kcoef<nkcoefs; kcoef++)
	fscanf(fp, "%le", &data->pespline[mbin].kspline[pcoef].coefs[kcoef]);
    splinealloc(&data->pespline[mbin].pspline, 1, npcoefs, 4,
		data->pespline[mbin].pmin, data->pespline[mbin].pmax);
    splinealloc(&data->pespline[mbin].ddkspline, 1, npcoefs, 4,
		data->pespline[mbin].pmin, data->pespline[mbin].pmax);
  } /* end loop mbin */

  /* Close input file */
  fclose(fp);

  return data;
}
#else
/******************************************************************************/
/* MPI version of spline coefficient reader: root reads, bcasts. */
spline3d *readspline(const char *fname)
{
  FILE     *fp;
  spline3d *data;
  double   emin;
  int      myrank, mbin, pcoef, kcoef, npcoefs, nkcoefs;

  /* Allocate storage */
  data = (spline3d *)malloc(sizeof(spline3d));
  data->norm = data->pmin = 1.0;  data->pmax = data->emax = -1.0;  data->cmin=1.0e+9;

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (!myrank) { /* Only root process reads file */

    /* Open input file */
    if ((fp = fopen(fname, "r")) == NULL) {
      fprintf(stderr, "Error: could not open file %s for reading.\n", fname);
      exit(1);
    }

    /* Read data */
    fscanf(fp, "%d", &data->nmubins);
    if (data->nmubins < 1) exit(1);
    if ((data->mubounds = (double *)malloc((data->nmubins+1) * sizeof(double))) == NULL) {
      fputs("Out of memory in readspline.\n", stderr);
      exit(2);
    }
    data->mubounds[0] = 0.0;
    if ((data->pespline = (spline2d *)malloc(data->nmubins * sizeof(spline2d))) == NULL) {
      fputs("Out of memory in readspline.\n", stderr);
      exit(2);
    }
    for (mbin=0; mbin<data->nmubins; mbin++) {
      fscanf(fp, "%le", &data->mubounds[mbin+1]);
      if (data->mubounds[mbin] >= data->mubounds[mbin+1]) exit(1);
      fscanf(fp, "%le\t%le\t%le\t%le",
	     &data->pespline[mbin].cmin, &data->pespline[mbin].emax,
	     &data->pespline[mbin].pmin, &data->pespline[mbin].pmax);
      if (data->pespline[mbin].pmin < data->pmin) data->pmin = data->pespline[mbin].pmin;
      if (data->pespline[mbin].pmax > data->pmax) data->pmax = data->pespline[mbin].pmax;
      if (data->pespline[mbin].cmin < data->cmin) data->cmin = data->pespline[mbin].cmin;
      if (data->pespline[mbin].emax > data->emax) data->emax = data->pespline[mbin].emax;
      if (data->pespline[mbin].pmin >= data->pespline[mbin].pmax) exit(1);
      emin = data->pespline[mbin].cmin * data->mubounds[mbin];
      fscanf(fp, "%d\t%d", &npcoefs, &nkcoefs);
      if ((npcoefs < 4) || (nkcoefs < 4)) exit(1);
      data->pespline[mbin].npcoefs = npcoefs;
      data->pespline[mbin].nkcoefs = nkcoefs;
      splinealloc(&data->pespline[mbin].kspline, npcoefs, nkcoefs, 4,
		  emin, data->pespline[mbin].emax);
      for (pcoef=0; pcoef<npcoefs; pcoef++)
	for (kcoef=0; kcoef<nkcoefs; kcoef++)
	  fscanf(fp, "%le", &data->pespline[mbin].kspline[pcoef].coefs[kcoef]);
      splinealloc(&data->pespline[mbin].pspline, 1, npcoefs, 4,
		  data->pespline[mbin].pmin, data->pespline[mbin].pmax);
      splinealloc(&data->pespline[mbin].ddkspline, 1, npcoefs, 4,
		  data->pespline[mbin].pmin, data->pespline[mbin].pmax);
    } /* end loop mbin */

    /* Close input file */
    fclose(fp);
  }

  MPI_Bcast(&data->nmubins, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (myrank) {
    if ((data->mubounds = (double *)malloc((data->nmubins+1) * sizeof(double))) == NULL) {
      fputs("Out of memory in readspline.\n", stderr);
      exit(2);
    }
    if ((data->pespline = (spline2d *)malloc(data->nmubins * sizeof(spline2d))) == NULL) {
      fputs("Out of memory in readspline.\n", stderr);
      exit(2);
    }
  }
  MPI_Bcast(data->mubounds, data->nmubins+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  for (mbin=0; mbin<data->nmubins; mbin++) {
    MPI_Bcast(&data->pespline[mbin].npcoefs, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&data->pespline[mbin].nkcoefs, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&data->pespline[mbin].cmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&data->pespline[mbin].emax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&data->pespline[mbin].pmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&data->pespline[mbin].pmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (myrank) {
      emin = data->pespline[mbin].cmin * data->mubounds[mbin];
      splinealloc(&data->pespline[mbin].kspline,
		  data->pespline[mbin].npcoefs, data->pespline[mbin].nkcoefs, 4,
		  emin, data->pespline[mbin].emax);
      splinealloc(&data->pespline[mbin].pspline,
		  1, data->pespline[mbin].npcoefs, 4,
		  data->pespline[mbin].pmin, data->pespline[mbin].pmax);
      splinealloc(&data->pespline[mbin].ddkspline,
		  1, data->pespline[mbin].npcoefs, 4,
		  data->pespline[mbin].pmin, data->pespline[mbin].pmax);
    }
    for (pcoef=0; pcoef<data->pespline[mbin].npcoefs; pcoef++)
      MPI_Bcast(data->pespline[mbin].kspline[pcoef].coefs,
		data->pespline[mbin].nkcoefs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  MPI_Bcast(&data->pmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&data->pmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&data->cmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&data->emax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  return data;
}
#endif

#ifdef PDF_NEAREST_NEIGHBOR
/******************************************************************************/
/* Note: uses nearest-neighbor interpolation in mu direction. */
double pdfnn(const double pphi, const double mu, const double ke, spline3d *data)
{
  double pdf;
  int    mbin, pcoef, pstart;

  /* Find the mu bin containing this point */
  if (mu < 0.0) return 0.0;
  for (mbin=0; mbin<data->nmubins; mbin++)
    if (mu < data->mubounds[mbin+1]) break;
  if (mbin == data->nmubins) return 0.0;

  /* Construct the vector of pphi coefficients here */
  pstart = getnzstart(data->pespline[mbin].pspline, pphi);
  for (pcoef=pstart; pcoef<pstart+4; pcoef++)
    getsplineval(&(data->pespline[mbin].kspline[pcoef]), ke,
		 &(data->pespline[mbin].pspline->coefs[pcoef]));

  /* Interpolate to get pdf val */
  getsplineval(data->pespline[mbin].pspline, pphi, &pdf);
  if (pdf < 0.0) return 0.0;
  if (data->logbins) pdf = exp(pdf);
  return data->norm*pdf;
}

/******************************************************************************/
/* Note: uses nearest-neighbor interpolation in mu direction. */
double pdfnnder(const double pphi, const double mu, const double ke, double *derivs,
		spline3d *data)
{
  double ys[2], factor=1.0;
  int    mbin, pcoef, pstart;

  /* Find the mu bin containing this point */
  if (mu < 0.0) return *derivs = derivs[1] = 0.0;
  for (mbin=0; mbin<data->nmubins; mbin++)
    if (mu < data->mubounds[mbin+1]) break;
  if (mbin == data->nmubins) return *derivs = derivs[1] = 0.0;

  /* Construct the vector of pphi coefficients here */
  pstart = getnzstart(data->pespline[mbin].pspline, pphi);
  for (pcoef=pstart; pcoef<pstart+4; pcoef++) {
    getsplinederiv(&(data->pespline[mbin].kspline[pcoef]), ke, ys);
    data->pespline[mbin].pspline->coefs[pcoef] = *ys;
    data->pespline[mbin].ddkspline->coefs[pcoef] = ys[1];
  }

  /* Interpolate to get pdf, deriv vals */
  getsplinederiv(data->pespline[mbin].pspline, pphi, ys);
  if (*ys < 0.0) return *derivs = derivs[1] = 0.0;
  if (data->logbins) factor = (ys[0] = exp(ys[0]));
  *derivs = data->norm*factor*ys[1];
  getsplineval(data->pespline[mbin].ddkspline, pphi, derivs+1);
  derivs[1] *= data->norm*factor;

  return data->norm*ys[0];
}

#else
/******************************************************************************/
/* Note: uses linear interpolation in mu direction. */
double pdflin(const double pphi, const double mu, const double ke, spline3d *data)
{
  double mean=0.0, oldmean=0.0, factor0, pdf0, pdf1;
  int    mbin, pcoef, pstart;

  /* Find the mu bin containing this point */
  for (mbin=0; mbin<data->nmubins; mbin++) {
    mean = 0.5*(data->mubounds[mbin] + data->mubounds[mbin+1]);
    if (mu < mean) break;
    oldmean = mean;
  }
  if (mbin == data->nmubins) { /* Last half-bin: interp to 0 at upper bndry */
    factor0 = (data->mubounds[mbin] - mu)/(data->mubounds[mbin] - mean);
    if (factor0 < 0.0) factor0 = 0.0;
    pdf1 = 0.0;
  } else { /* Not last bin */
    /* Construct the vector of pphi coefficients here */
    pstart = getnzstart(data->pespline[mbin].pspline, pphi);
    for (pcoef=pstart; pcoef<pstart+4; pcoef++)
      getsplineval(&(data->pespline[mbin].kspline[pcoef]), ke,
		   &(data->pespline[mbin].pspline->coefs[pcoef]));

    /* Interpolate to get pdf val */
    getsplineval(data->pespline[mbin].pspline, pphi, &pdf1);

    if (!mbin) { /* 1st half-bin: do not interpolate */
      if (pdf1 < 0.0) return 0.0;
      return data->norm*pdf1;
    }

    factor0 = (mean - mu)/(mean - oldmean);
  }

  /* Construct the vector of pphi coefficients here */
  pstart = getnzstart(data->pespline[mbin-1].pspline, pphi);
  for (pcoef=pstart; pcoef<pstart+4; pcoef++)
    getsplineval(&(data->pespline[mbin-1].kspline[pcoef]), ke,
		 &(data->pespline[mbin-1].pspline->coefs[pcoef]));

  /* Interpolate to get pdf val */
  getsplineval(data->pespline[mbin-1].pspline, pphi, &pdf0);

  pdf1 = factor0*pdf0 + (1.0 - factor0)*pdf1;
  if (pdf1 < 0.0) return 0.0;
  return data->norm*pdf1;
}

/******************************************************************************/
/* Note: uses linear interpolation in mu direction. */
double pdflinder(const double pphi, const double mu, const double ke, double *derivs,
		 spline3d *data)
{
  double ys[2];
  double mean=0.0, oldmean=0.0, factor0, pdf0, ddk0, ddp0, pdf1;
  int    mbin, pcoef, pstart;

  /* Find the mu bin containing this point */
  for (mbin=0; mbin<data->nmubins; mbin++) {
    mean = 0.5*(data->mubounds[mbin] + data->mubounds[mbin+1]);
    if (mu < mean) break;
    oldmean = mean;
  }
  if (mbin == data->nmubins) {
    factor0 = (data->mubounds[mbin] - mu)/(data->mubounds[mbin] - mean);
    if (factor0 < 0.0) factor0 = 0.0;
    pdf1 = 0.0;
  } else {
    factor0 = (mean - mu)/(mean - oldmean);

    /* Construct the vector of pphi, pphi-prime coefficients here */
    pstart = getnzstart(data->pespline[mbin].pspline, pphi);
    for (pcoef=pstart; pcoef<pstart+4; pcoef++) {
      getsplinederiv(&(data->pespline[mbin].kspline[pcoef]), ke, ys);
      data->pespline[mbin].pspline->coefs[pcoef] = *ys;
      data->pespline[mbin].ddkspline->coefs[pcoef] = ys[1];
    }

    /* Interpolate to get pdf, deriv vals */
    getsplinederiv(data->pespline[mbin].pspline, pphi, ys);
    pdf1 = *ys;  *derivs = ys[1];
    getsplineval(data->pespline[mbin].ddkspline, pphi, derivs+1);

    if (!mbin) {
      if (pdf1 < 0.0) return *derivs = derivs[1] = 0.0;
      derivs[0] *= data->norm;  derivs[1] *= data->norm;
      return data->norm*pdf1;
    }
  }

  /* Construct the vector of pphi, pphi-prime coefficients here */
  pstart = getnzstart(data->pespline[mbin-1].pspline, pphi);
  for (pcoef=pstart; pcoef<pstart+4; pcoef++) {
    getsplinederiv(&(data->pespline[mbin-1].kspline[pcoef]), ke, ys);
    data->pespline[mbin-1].pspline->coefs[pcoef] = *ys;
    data->pespline[mbin-1].ddkspline->coefs[pcoef] = ys[1];
  }

  /* Interpolate to get pdf, deriv vals */
  getsplinederiv(data->pespline[mbin-1].pspline, pphi, ys);
  pdf0 = *ys;  ddp0 = ys[1];
  getsplineval(data->pespline[mbin-1].ddkspline, pphi, &ddk0);

  *derivs = factor0*ddp0 + (1.0 - factor0)*(*derivs);
  derivs[1] = factor0*ddk0 + (1.0 - factor0)*derivs[1];
  pdf1 = factor0*pdf0 + (1.0 - factor0)*pdf1;
  if (pdf1 < 0.0) return *derivs = derivs[1] = 0.0;
  derivs[0] *= data->norm;  derivs[1] *= data->norm;
  return data->norm*pdf1;
}
#endif

/******************************************************************************/
double integratedV(spline3d *data,
		   const double pi, const double pf, const int np,
		   const double mi, const double mf, const int nm,
		   const double ki, const double kf, const int nk)
{
  double mu = mi, pp, ke;
  double dm = (mf - mi)/(nm - 1);
  double dp = (pf - pi)/(np - 1);
  double dk = (kf - ki)/(nk - 1);
  double sum=0.0;
  int    ip,im,ik;

  for (im=0; im<nm; im++, mu+=dm)
    for (ip=0, pp=pi; ip<np; ip++, pp+=dp)
      for (ik=0, ke=ki; ik<nk; ik++, ke+=dk)
	sum += PDFFN(pp, mu, ke, data);

  return dm*dp*dk*sum;
}

/******************************************************************************/
void freespline2d(spline2d *thespline)
{
  int i;

  splinedelete(thespline->pspline);
  free(thespline->pspline);
  splinedelete(thespline->ddkspline);
  free(thespline->ddkspline);
  for (i=0; i<thespline->npcoefs; i++)
    splinedelete(thespline->kspline + i);
  free(thespline->kspline);
}

/******************************************************************************/
void freespline3d(spline3d *thespline)
{
  int i;

  free(thespline->mubounds);
  for (i=0; i<thespline->nmubins; i++)
    freespline2d(thespline->pespline + i);
  free(thespline->pespline);
  free(thespline);
}
