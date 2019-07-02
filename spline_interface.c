/* module load gsl; gcc -O -c spline_interface.c */
#include "spline_interface.h"

/******************************************************************************/
/* Create a 1D spline to fit a given 1D dataset with specified bounds, using  */
/* GNU scientific library (gsl) routines.                                     */
spline *createBspline1d(const double *xarr, const double *yarr, const int ndata,
			const int ncoefs, const int korder,
			const double lbound, const double ubound,
			const unsigned long flags)
{
  spline                        *thespline;
  gsl_multifit_linear_workspace *fwork;
  gsl_matrix                    *fitmat, *cov;
  gsl_vector                    *coefs, *yvec;
  double                         chisq;
  int                            jj, icol;

  /* Error check */
  if (korder < 1) {
    fprintf(stderr, "Spline order %d too low in createBspline1d.\n", korder);
    return NULL;
  }
  if (ncoefs >= ndata) {
    fputs("ncoefs >= ndata; fit is underdetermined.\n", stderr);
    return NULL;
  }

  /* Determine bounds of dataset */
  thespline = (spline *)malloc(sizeof(spline));
  thespline->k = korder;
  if (flags & LOWER_BOUND_ZERO)
    thespline->lbound = 0.0;
  else
    thespline->lbound = lbound;
  thespline->ubound = ubound;

  /* Set up uniform knots */
  thespline->bwork = gsl_bspline_alloc(korder, ncoefs + 2 - korder);
  gsl_bspline_knots_uniform(thespline->lbound, thespline->ubound, thespline->bwork);
  thespline->dwork = gsl_bspline_deriv_alloc(korder);
  thespline->bval = gsl_vector_alloc(korder);
  thespline->nzwork = gsl_matrix_alloc(korder, 2);

  /* Construct matrix, LHS for linear least-squares fit */
  yvec = gsl_vector_alloc(ndata);
  coefs = gsl_vector_alloc(ncoefs);
  fitmat = gsl_matrix_alloc(ndata, ncoefs);
  for (jj=0; jj<ndata; jj++) {
    /* Insert y data point into GSL vector */
    gsl_vector_set(yvec, jj, yarr[jj]);

    /* Evaluate B-spline basis functions at this point */
    gsl_bspline_eval(xarr[jj], coefs, thespline->bwork);

    /* Insert them into the appropriate row of the fitting matrix */
    for (icol=0; icol<ncoefs; icol++)
      gsl_matrix_set(fitmat, jj, icol, gsl_vector_get(coefs, icol));
  } /* end loop jj */

  /* Perform linear least-squares fit */
  fwork = gsl_multifit_linear_alloc(ndata, ncoefs);
  cov = gsl_matrix_alloc(ncoefs, ncoefs);
  gsl_multifit_linear(fitmat, yvec, coefs, cov, &chisq, fwork);
  gsl_matrix_free(fitmat);  gsl_matrix_free(cov);
  gsl_vector_free(yvec);  gsl_multifit_linear_free(fwork);

  /* Evaluate goodness of fit */
  chisq = chisq/(ndata - ncoefs);
  /* fprintf(stderr, "chi-squared/dof = %le\n", chisq);
     if (chisq > 1.0e-4) fputs("Warning: spline fit is poor.\n", stderr); */

  /* Copy spline coefficients to permanent data structure */
  thespline->ncoefs = ncoefs;
  if ((thespline->coefs = (double *)malloc(ncoefs * sizeof(double))) == NULL) {
    fputs("Insufficient memory for spline coefficients in createBspline1d.\n", stderr);
    free(thespline);
    return NULL;
  }
  for (icol=0; icol<ncoefs; icol++)
    thespline->coefs[icol] = gsl_vector_get(coefs, icol);
  gsl_vector_free(coefs);

  return thespline;
}

/******************************************************************************/
/* Allocate space for an array of 1D splines of a specified size and bounds;  */
/* initialize knots but not coefficients.                                     */
void splinealloc(spline **splinearr, const int nsplines, const int ncoefs,
		 const int korder, const double lbound, const double ubound)
{
  int ispline;

  if ((*splinearr = (spline *)malloc(nsplines * sizeof(spline))) == NULL) {
    fputs("Out of memory in splinealloc.\n", stderr);
    return ;
  }

  for (ispline=0; ispline<nsplines; ispline++) {
    (*splinearr)[ispline].k = korder;
    (*splinearr)[ispline].ncoefs = ncoefs;
    (*splinearr)[ispline].bwork = gsl_bspline_alloc(korder, ncoefs + 2 - korder);
    (*splinearr)[ispline].dwork = gsl_bspline_deriv_alloc(korder);
    (*splinearr)[ispline].bval = gsl_vector_alloc(korder);
    (*splinearr)[ispline].nzwork = gsl_matrix_alloc(korder, 2);

    (*splinearr)[ispline].coefs = (double *)malloc(ncoefs * sizeof(double));
    if ((*splinearr)[ispline].coefs == NULL) {
      fputs("Out of memory in splinealloc.\n", stderr);
      return ;
    }

    /* Set up bounds, knots */
    (*splinearr)[ispline].lbound = lbound;
    (*splinearr)[ispline].ubound = ubound;
    gsl_bspline_knots_uniform(lbound, ubound, (*splinearr)[ispline].bwork);
  } /* end loop ispline */
}

/******************************************************************************/
/* Delete all internal allocated storage for a single spline set up with      */
/* of the two above routines, createBspline1d() or splinealloc(). This should */
/* be called before the spline itself is freed with free() to avoid memory    */
/* leaks. */
void splinedelete(spline *thespline)
{
  gsl_bspline_free(thespline->bwork);
  gsl_bspline_deriv_free(thespline->dwork);
  gsl_vector_free(thespline->bval);
  gsl_matrix_free(thespline->nzwork);
  free(thespline->coefs);
}

/******************************************************************************/
/* Use the given spline to interpolate the function y(x).                     */
void getsplineval(spline *thespline, const double x, double *y)
{
  size_t istart, iend;
  int    i;

  *y = 0.0;
  if ((x < thespline->lbound) || (x > thespline->ubound)) return;

  gsl_bspline_eval_nonzero(x, thespline->bval, &istart, &iend, thespline->bwork);
 
  for (i=0; i<thespline->k; i++) 
    *y += thespline->coefs[istart+i] * gsl_vector_get(thespline->bval,i);
}

/******************************************************************************/
/* Function and its 1st derivative at x */
void getsplinederiv(spline *thespline, const double x, double *ys)
{
  double c;
  size_t istart, iend;
  int    i;

  *ys = ys[1] = 0.0;
  if ((x < thespline->lbound) || (x > thespline->ubound)) return;

  gsl_bspline_deriv_eval_nonzero(x, 1, thespline->nzwork, &istart, &iend,
				 thespline->bwork, thespline->dwork);

  for (i=0; i<thespline->k; i++) {
    c = thespline->coefs[istart+i];
    ys[0] += c * gsl_matrix_get(thespline->nzwork,i,0);
    ys[1] += c * gsl_matrix_get(thespline->nzwork,i,1);
  }
}

/******************************************************************************/
/* Determine the index of the first basis spline that will have a nonzero     */
/* value when computing the function or its nonzero values at x.              */
int getnzstart(spline *s, const double x)
{
  size_t i, f;
  if (x < s->lbound) return 0;
  gsl_bspline_eval_nonzero(x, s->bval, &i, &f, s->bwork);
  return (int)i;
}
