#ifndef PARTICLESPLINES_H
#define PARTICLESPLINES_H

#include "spline_interface.h"

/* Type definitions */
typedef struct {
  spline *kspline, *pspline, *ddkspline;
  double  cmin, emax, pmin, pmax;
  int     npcoefs, nkcoefs;
} spline2d;

typedef struct {
  spline2d *pespline;
  double   *mubounds;
  double    norm, pmin, pmax, cmin, emax;
  int       logbins, nmubins, iorder;
} spline3d;

/* Function prototypes */
spline3d *readspline(const char *);
inline int mubin(spline3d *, const double);
double pdf2d(spline2d *, const double, const double, const int);
double pdf2dder(spline2d *, const double, const double, const int, double *);
double pdfnn(const double, const double, const double, spline3d *);
double pdfnnder(const double, const double, const double, double *,
		spline3d *);
double pdflin(const double, const double, const double, spline3d *);
double pdflinder(const double, const double, const double, double *,
		 spline3d *);
double integratedV(spline3d *, const double, const double, const int,
		   const double, const double, const int,
		   const double, const double, const int);
void freespline2d(spline2d *);
void freespline3d(spline3d *);

#endif
