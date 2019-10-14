#ifndef SPLINE_INTERFACE_H
#define SPLINE_INTERFACE_H

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>

/* Macros */
#define LOWER_BOUND_ZERO   1

/* Data type definition */
typedef struct {
  gsl_bspline_workspace       *bwork; /* gsl internal spline workspace, contains knot info */
  gsl_vector                  *bval; /* Workspace for basis spline evaluation at x */
  gsl_matrix                  *nzwork; /* Workspace for basis spline, deriv evaluation */
  double                      *coefs; /* Array of coefficients of basis splines */
  double                       lbound, ubound; /* Lower and upper bounds of splined domain */
  int                          k, ncoefs; /* Order of the spline+1 (4=cubic); # of basis pts */
} spline;

/* Function prototypes */
spline *createBspline1d(const double *xarr, const double *yarr, const int ndata,
			const int ncoefs, const int korder,
			const double, const double,
			const unsigned long flags);
void splinealloc(spline **sa, const int ns, const int nc, const int k,
		 const double lb, const double ub);
void splinedelete(spline *);
void getsplineval(spline *thespline, const double x, double *y);
void getsplinederiv(spline *thespline, const double x, double *y);
int  getnzstart(spline *, const double);

#endif
