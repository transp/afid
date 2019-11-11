#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "particleSplines.h"
#ifdef PARTICLE_MPI
#include "mpi.h"
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

  fputs("Setting up splines for interpolation.\n", stderr);

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

////////////////////////////////////////////////////////////////////////////////
int pspline_set_order_(int *iorder)
{
  nsplinedata->iorder = *iorder;
  psplinedata->iorder = *iorder;

  switch(*iorder) {
  case 0:
    fputs("Setting interpolation type to nearest-neighbor.\n", stderr);
    break;
  case 1:
    fputs("Setting interpolation type to linear.\n", stderr);
    break;
  default:
    fprintf(stderr, "Warning: unrecognized interpolation type %d requested.\n", *iorder);
    return 1;
  }

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
  spline3d *data;

  data = (*sgnv > 0) ? psplinedata : nsplinedata;

  switch(data->iorder) {
  case 0:
    *f = pdfnn(*pphi, *mu, *ke, data);
    break;
  case 1:
    *f = pdflin(*pphi, *mu, *ke, data);
    break;
  default:
    fprintf(stderr, "Unsupported interpolation order %d in getpdf\n", data->iorder);
    *f = 0.0;
    return 1;
  }

  return 0;
}

/******************************************************************************/
/* Evaluate f(p_phi, mu, E) and its 1st derivatives with respect to P and E */
/*  by 2D cubic B-spline interpolation in P,E */
int getpdfd_(double *pphi, double *mu, double *ke, int *sgnv, double *f,
	    double *dfdp, double *dfdE)
{
  spline3d *data;
  double drv[2];

  data = (*sgnv > 0) ? psplinedata : nsplinedata;

  switch(data->iorder) {
  case 0:
    *f = pdfnnder(*pphi, *mu, *ke, drv, data);
    break;
  case 1:
    *f = pdflinder(*pphi, *mu, *ke, drv, data);
    break;
  default:
    fprintf(stderr, "Unsupported interpolation order %d in getpdfd\n", data->iorder);
    *f = *dfdp = *dfdE = 0.0;
    return 1;
  }

  *dfdp = drv[0];  *dfdE = drv[1];

  return 0;
}

/******************************************************************************/
/* Evaluate f(p_phi, mu, E) and its 1st derivative with respect to E          */
/*  by 2D cubic B-spline interpolation in P,E */
int getdfde_(double *pphi, double *mu, double *ke, int *sgnv, double *f,
	     double *dfdE)
{
  spline3d *data;

  data = (*sgnv > 0) ? psplinedata : nsplinedata;

  switch(data->iorder) {
  case 0:
    *f = pdfnndE(*pphi, *mu, *ke, dfdE, data);
    break;
    //case 1:
    //*f = pdflinder(*pphi, *mu, *ke, drv, data);
    //break;
  default:
    fprintf(stderr, "Unsupported interpolation order %d in getdfdE\n", data->iorder);
    *f = *dfdE = 0.0;
    return 1;
  }

  return 0;
}

/******************************************************************************/
/* Evaluate f(p_phi, mu, E) and its 1st derivative with respect to P_phi      */
/*  by 2D cubic B-spline interpolation in P,E */
int getdfdp_(double *pphi, double *mu, double *ke, int *sgnv, double *f,
	     double *dfdP)
{
  spline3d *data;

  data = (*sgnv > 0) ? psplinedata : nsplinedata;

  switch(data->iorder) {
  case 0:
    *f = pdfnndP(*pphi, *mu, *ke, dfdP, data);
    break;
    //case 1:
    //*f = pdflinder(*pphi, *mu, *ke, drv, data);
    //break;
  default:
    fprintf(stderr, "Unsupported interpolation order %d in getdfdE\n", data->iorder);
    *f = *dfdP = 0.0;
    return 1;
  }

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
  data->logbins = 0;  data->iorder = 0;

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

/******************************************************************************/
// Find the mu bin containing the specified mu value.
inline int mubin(spline3d *data, const double mu)
{
  int mbin;

  for (mbin=-1; mbin<data->nmubins; mbin++)
    if (mu < data->mubounds[mbin+1]) break;

  return mbin;
}

/******************************************************************************/
// Evaluate the 2D spline at the specified pphi, ke point.
double pdf2d(spline2d *data, const double pphi, const double ke, const int logbins)
{
  double pdf;
  int    pcoef, pstart;

  // Construct the vector of pphi coefficients for f at the specified ke.
  pstart = getnzstart(data->pspline, pphi);
  for (pcoef=pstart; pcoef<pstart+4; pcoef++)
    getsplineval(data->kspline + pcoef, ke, &(data->pspline->coefs[pcoef]));

  // Interpolate to get pdf value.
  getsplineval(data->pspline, pphi, &pdf);
  if (logbins) pdf = exp(pdf);
  else if (pdf < 0.0) return 0.0; // Rectify negative f values
  return pdf;
}

/******************************************************************************/
// Evaluate the 2D spline and both its first partial derivatives at
//  the specified pphi, ke point.
double pdf2dder(spline2d *data, const double pphi, const double ke,
		const int logbins, double *fprime)
{
  double ys[2], factor=1.0;
  int    pcoef, pstart;

  // Construct the vector of f, df/dk pphi coefficients here
  pstart = getnzstart(data->pspline, pphi);
  for (pcoef=pstart; pcoef<pstart+4; pcoef++) {
    getsplinederiv(data->kspline + pcoef, ke, ys);
    data->pspline->coefs[pcoef]   = ys[0];
    data->ddkspline->coefs[pcoef] = ys[1];
  }

  // Interpolate to get pdf, derivative values
  //  If logbins then g == ln(f) else g == f
  getsplinederiv(data->pspline, pphi, ys); // ys[0] = g;  ys[1] = dg/dpphi
  if (logbins) factor = (ys[0] = exp(ys[0])); // f
  else if (ys[0] < 0.0) return fprime[0] = fprime[1] = 0.0; // Rectify f < 0
  fprime[0] = factor*ys[1]; // df/pphi
  getsplineval(data->ddkspline, pphi, fprime+1);
  fprime[1] *= factor;      // df/dke

  return ys[0];
}

/******************************************************************************/
// Evaluate the 2D spline and its first partial derivative with respect to E
//  at the specified pphi, ke point.
double pdf2ddfdE(spline2d *data, const double pphi, const double ke,
		 const int logbins, double *dfdE)
{
  double ys[2], factor=1.0;
  int    pcoef, pstart;

  // Construct the vector of f, df/dk pphi coefficients here
  pstart = getnzstart(data->pspline, pphi);
  for (pcoef=pstart; pcoef<pstart+4; pcoef++) {
    getsplinederiv(data->kspline + pcoef, ke, ys);
    data->pspline->coefs[pcoef]   = ys[0];
    data->ddkspline->coefs[pcoef] = ys[1];
  }

  // Interpolate to get pdf, derivative values
  //  If logbins then g == ln(f) else g == f
  getsplineval(data->pspline, pphi, ys); // ys[0] = g
  if (logbins) factor = (ys[0] = exp(ys[0])); // f
  else if (ys[0] < 0.0) return *dfdE = 0.0; // Rectify f < 0
  getsplineval(data->ddkspline, pphi, dfdE);
  *dfdE *= factor;      // df/dke

  return ys[0];
}

/******************************************************************************/
// Evaluate the 2D spline and its first partial derivative with respect to P_phi
//  at the specified pphi, ke point.
double pdf2ddfdP(spline2d *data, const double pphi, const double ke,
		 const int logbins, double *dfdP)
{
  double ys[2], factor=1.0;
  int    pcoef, pstart;

  // Construct the vector of f, df/dk pphi coefficients here
  pstart = getnzstart(data->pspline, pphi);
  for (pcoef=pstart; pcoef<pstart+4; pcoef++)
    getsplineval(data->kspline + pcoef, ke, &(data->pspline->coefs[pcoef]));

  // Interpolate to get pdf, derivative values
  //  If logbins then g == ln(f) else g == f
  getsplinederiv(data->pspline, pphi, ys); // ys[0] = g;  ys[1] = dg/dpphi
  if (logbins) factor = (ys[0] = exp(ys[0])); // f
  else if (ys[0] < 0.0) return *dfdP = 0.0; // Rectify f < 0
  *dfdP = factor*ys[1]; // df/pphi

  return ys[0];
}

/******************************************************************************/
// Evaluate f(P,mu,E) using nearest-neighbor interpolation in mu
double pdfnn(const double pphi, const double mu, const double ke, spline3d *data)
{
  int mbin;

  // Find the mu bin containing this point
  mbin = mubin(data, mu);
  if ((mbin < 0) || (mbin >= data->nmubins)) return 0.0; // out of range!

  // Return the normalized spline value
  return data->norm*pdf2d(data->pespline + mbin, pphi, ke, data->logbins);
}

/******************************************************************************/
// Evaluate f(P,mu,E), df/dP, df/dE using nearest-neighbor interpolation in mu
double pdfnnder(const double pphi, const double mu, const double ke, double *derivs,
		spline3d *data)
{
  double fn;
  int    mbin;

  // Find the mu bin containing this point
  mbin = mubin(data, mu);
  if ((mbin < 0) || (mbin >= data->nmubins)) // out of range!
    return derivs[0] = derivs[1] = 0.0;

  // Return the normalized spline, derivative values
  fn = pdf2dder(data->pespline + mbin, pphi, ke, data->logbins, derivs);
  derivs[0] *= data->norm;  derivs[1] *= data->norm;
  return data->norm*fn;
}

/******************************************************************************/
// Evaluate f(P,mu,E), df/dE using nearest-neighbor interpolation in mu
double pdfnndE(const double pphi, const double mu, const double ke, double *dE,
	       spline3d *data)
{
  double fn;
  int    mbin;

  // Find the mu bin containing this point
  mbin = mubin(data, mu);
  if ((mbin < 0) || (mbin >= data->nmubins)) // out of range!
    return *dE = 0.0;

  // Return the normalized spline, derivative values
  fn = pdf2ddfdE(data->pespline + mbin, pphi, ke, data->logbins, dE);
  *dE *= data->norm;
  return data->norm*fn;
}

/******************************************************************************/
// Evaluate f(P,mu,E), df/dPphi using nearest-neighbor interpolation in mu
double pdfnndP(const double pphi, const double mu, const double ke, double *dP,
	       spline3d *data)
{
  double fn;
  int    mbin;

  // Find the mu bin containing this point
  mbin = mubin(data, mu);
  if ((mbin < 0) || (mbin >= data->nmubins)) // out of range!
    return *dP = 0.0;

  // Return the normalized spline, derivative values
  fn = pdf2ddfdP(data->pespline + mbin, pphi, ke, data->logbins, dP);
  *dP *= data->norm;
  return data->norm*fn;
}

/******************************************************************************/
// Evaluate f(P,mu,E) using linear interpolation in mu
double pdflin(const double pphi, const double mu, const double ke, spline3d *data)
{
  double midpt, factor;
  int    mbin, endflg=0;

  // Find the mu bin containing this point
  mbin = mubin(data, mu);
  if ((mbin < 0) || (mbin >= data->nmubins)) return 0.0; // out of range!

  // Determine which two bins will be averaged, interpolation factor
  midpt = 0.5*(data->mubounds[mbin] + data->mubounds[mbin+1]);
  factor = 1.0 - 0.5*(mu - midpt)/(data->mubounds[mbin+1] - midpt);
  if (factor > 1.0) { factor -= 1.0;  mbin--;}

  // Special cases at ends
  if (mbin < 0) { endflg=1;  mbin=0;}
  else if (mbin >= data->nmubins-1) { endflg=1;  mbin = data->nmubins-1;}
  if (endflg) return data->norm*pdf2d(data->pespline + mbin, pphi, ke, data->logbins);

  // General case
  return data->norm*(factor*pdf2d(data->pespline + mbin, pphi, ke, data->logbins) +
		     (1.0-factor)*pdf2d(data->pespline + mbin+1, pphi, ke, data->logbins));
}

/******************************************************************************/
// Evaluate f(P,mu,E), df/dP, df/dE using linear interpolation in mu
double pdflinder(const double pphi, const double mu, const double ke, double *derivs,
		 spline3d *data)
{
  double midpt, factor, fn1, fn2, tmp[2];
  int    mbin, endflg=0;

  // Find the mu bin containing this point
  mbin = mubin(data, mu);
  if ((mbin < 0) || (mbin >= data->nmubins)) // out of range!
    return derivs[0] = derivs[1] = 0.0;

  // Determine which two bins will be averaged, interpolation factor
  midpt = 0.5*(data->mubounds[mbin] + data->mubounds[mbin+1]);
  factor = 1.0 - 0.5*(mu - midpt)/(data->mubounds[mbin+1] - midpt);
  if (factor > 1.0) { factor -= 1.0;  mbin--;}

  // Special cases at ends
  if (mbin < 0) { endflg=1;  mbin=0;}
  else if (mbin >= data->nmubins-1) { endflg=1;  mbin = data->nmubins-1;}
  if (endflg) {
    fn1 = pdf2dder(data->pespline + mbin, pphi, ke, data->logbins, derivs);
    derivs[0] *= data->norm;  derivs[1] *= data->norm;
    return data->norm*fn1;
  }

  // General case
  fn1 = pdf2dder(data->pespline + mbin,   pphi, ke, data->logbins, derivs);
  fn2 = pdf2dder(data->pespline + mbin+1, pphi, ke, data->logbins, tmp);
  derivs[0] = data->norm*(factor*derivs[0] + (1.0-factor)*tmp[0]);
  derivs[1] = data->norm*(factor*derivs[1] + (1.0-factor)*tmp[1]);
  return data->norm*(factor*fn1 + (1.0-factor)*fn2);
}

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
	sum += pdfnn(pp, mu, ke, data);

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
