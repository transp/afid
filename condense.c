#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <dirent.h>
#include <math.h>
#include "netcdf.h"

/* Type definitions */
typedef struct {
  double *mass, *charge;
  int     nnbi;
} species_info;

/* Global variable(s) */
char ppat[256]; /* Search pattern for particle data file names */

/* Function prototype(s) */
int filter(const struct dirent *);
int get_species_info(const char *, species_info *);

/******************************************************************************/
void usage(const char *exname)
{
  fprintf(stderr, "%s: condense NUBEAM particle output for binning.\n",
	  exname);
  fprintf(stderr, "\nUsage:\n%s <file identifier> [output directory path]\n\n", exname);
  exit(0);
}

/******************************************************************************/
int main(int argc, char *argv[])
{
  char            path[256];
  struct dirent **namelist;
  double         *lambda, *Rcoord, *zcoord, *muval, *pphi, *modv, *wt;
  double          Rmin=1e+30, Rmax=-1e+30, zmin=1e+30, zmax=-1e+30, vmax=-1e+30;
  double          minwt=1.0e+80, maxwt=0.0, vl, vr;
  nc_type         xtype;
  species_info    specinf;
  size_t          len0, len1, start[2], count[2], *last;
  int             dimids[2];
  int             ifile, ipart, ispec, nument, nptcls, ierr;
  int             ncid, ndims, nvars, ngatts, unlimdimid;
  int             mbid, lmbid, Rid, zid, specid, ptcid, muid, ppid, vid, wid;
  int             ncido, ptcdim, specdim, ppido, muido, Eido, wtido, lido, mido, cido;
  int             nnzf, nnzs, nnzt=0, first=1;
  const char      pitchname[] = "xksidy", muname[]="xmuay", ppname[]="pphiay";
  const char      vname[]="vay", Rname[]="rmjionay", zname[]="xzionay";
  const char      wname[] = "wghtay";
  const char      punits[]="kg m^2/s", munits[]="A m^2";
  signed char     cp;

  if (argc == 1) usage(*argv);

  /* Find files in current directory matching supplied pattern */
  sprintf(ppat, "%s_debug_nbi_ptcl_state_cpu", argv[1]);
  nument = scandir(".", &namelist, filter, alphasort);
  if (nument > 0)
    fprintf(stderr, "%d particle files found.\n", nument);
  else {
    fprintf(stderr, "No particle files found matching pattern %s.\n", argv[1]);
    return 0;
  }

  if (get_species_info(argv[1], &specinf)) {
    fputs("Assuming beam species is deuterium.\n", stderr);
    specinf.nnbi = 1;
    specinf.mass = (double *)malloc(sizeof(double));
    specinf.mass[0] = 3.343583e-27; /* deuteron mass */
    specinf.charge = (double *)malloc(sizeof(double));
    specinf.charge[0] = 1.6021765e-19; /* deuteron charge */
  }
  last = (size_t *)malloc(specinf.nnbi * sizeof(size_t));

  /* Create, initialize new NetCDF file to store condensed data */
  if (argc > 2)
    sprintf(path, "%s/%s_condensed.cdf", argv[2], argv[1]);
  else
    sprintf(path, "%s_condensed.cdf", argv[1]);
  ierr = nc_create(path, NC_CLOBBER, &ncido);
  if (ierr != NC_NOERR) {
    fprintf(stderr, "Error creating condensed file: %s.\n", nc_strerror(ierr));
    return 1;
  }
  ierr = nc_def_dim(ncido, "nptcl", NC_UNLIMITED, &ptcdim);
  ierr = nc_def_dim(ncido, "nspec", (size_t)specinf.nnbi, &specdim);
  dimids[0] = ptcdim;  dimids[1] = specdim;
  ierr = nc_def_var(ncido, "mass", NC_DOUBLE, 1, dimids+1, &mido);
  ierr = nc_put_att_text(ncido, mido, "units", 2, "kg");
  ierr = nc_def_var(ncido, "charge", NC_DOUBLE, 1, dimids+1, &cido);
  ierr = nc_put_att_text(ncido, cido, "units", 1, "C");
  ierr = nc_def_var(ncido, "pphi", NC_DOUBLE, 2, dimids, &ppido);
  ierr = nc_put_att_text(ncido, ppido, "units", strlen(punits), punits);
  ierr = nc_def_var(ncido, "mu", NC_DOUBLE, 2, dimids, &muido);
  ierr = nc_put_att_text(ncido, muido, "units", strlen(munits), munits);
  ierr = nc_def_var(ncido, "E", NC_DOUBLE, 2, dimids, &Eido);
  ierr = nc_put_att_text(ncido, Eido, "units", 1, "J");
  ierr = nc_def_var(ncido, "weight", NC_DOUBLE, 2, dimids, &wtido);
  ierr = nc_def_var(ncido, "sgn_v", NC_BYTE, 2, dimids, &lido);
  ierr = nc_enddef(ncido); /* leave define mode */

  /* Write species data */
  for (ispec=0; ispec<specinf.nnbi; ispec++) {
    *start = ispec;
    ierr = nc_put_var1_double(ncido, mido, start, &specinf.mass[ispec]);
    ierr = nc_put_var1_double(ncido, cido, start, &specinf.charge[ispec]);
    last[ispec] = 0L;
   }

  /* Loop over particle data files */
  for (ifile=0; ifile<nument; ifile++) {
    fprintf(stderr, "\nFile %d: %s\n", ifile, namelist[ifile]->d_name);

    /* Open the file */
    ierr = nc_open(namelist[ifile]->d_name, NC_NOWRITE, &ncid);
    if (ierr != NC_NOERR) {
      fprintf(stderr, "Error opening particle file: %s.\n", nc_strerror(ierr));
      continue;
    }

    /* Get overview of file contents */
    ierr = nc_inq(ncid, &ndims, &nvars, &ngatts, &unlimdimid);
    fprintf(stderr, "%d dimensions, %d variables, %d global attributes.\n",
	    ndims, nvars, ngatts);
    if (unlimdimid == -1) fputs("No unlimited dimension.\n", stderr);

    /* Check for presence of required variables */
    ierr = nc_inq_varid(ncid, "minb", &mbid);
    if (ierr != NC_NOERR) {
      fprintf(stderr, "minb: %s.\n", nc_strerror(ierr)); goto l10;
    }
    ierr = nc_inq_varid(ncid, pitchname, &lmbid);
    if (ierr != NC_NOERR) {
      fprintf(stderr, "%s: %s.\n", pitchname, nc_strerror(ierr)); goto l10;
    }
    ierr = nc_inq_varid(ncid, Rname, &Rid);
    if (ierr != NC_NOERR) {
      fprintf(stderr, "%s: %s.\n", Rname, nc_strerror(ierr)); goto l10;
    }
    ierr = nc_inq_varid(ncid, zname, &zid);
    if (ierr != NC_NOERR) {
      fprintf(stderr, "%s: %s.\n", zname, nc_strerror(ierr)); goto l10;
    }
    ierr = nc_inq_varid(ncid, muname, &muid);
    if (ierr != NC_NOERR) {
      fprintf(stderr, "%s: %s.\n", muname, nc_strerror(ierr)); goto l10;
    }
    ierr = nc_inq_varid(ncid, ppname, &ppid);
    if (ierr != NC_NOERR) {
      fprintf(stderr, "%s: %s.\n", ppname, nc_strerror(ierr)); goto l10;
    }
    ierr = nc_inq_varid(ncid, vname, &vid);
    if (ierr != NC_NOERR) {
      fprintf(stderr, "%s: %s.\n", vname, nc_strerror(ierr)); goto l10;
    }
    ierr = nc_inq_varid(ncid, wname, &wid);
    if (ierr != NC_NOERR) {
      fprintf(stderr, "%s: %s.\n", wname, nc_strerror(ierr)); goto l10;
    }

    /* Get variable dimensions */
    ierr = nc_get_var1_int(ncid, mbid, 0, &nptcls);
    if (ierr != NC_NOERR) {
      fprintf(stderr, "Error: %s.\n", nc_strerror(ierr)); goto l10;
    }
    fprintf(stderr, "%d particles in file.\n", nptcls);
    ierr = nc_inq_vartype(ncid, lmbid, &xtype);
    fprintf(stderr, "Particle data type is %d.\n", (int)xtype);
    if (xtype == NC_DOUBLE) fputs("(double precision)\n", stderr);
    else {
      if (xtype == NC_FLOAT) fputs("(single precision)\n", stderr);
      goto l10;
    }
    ierr = nc_inq_varndims(ncid, lmbid, &ndims);
    fprintf(stderr, "Particle data has %d dimension(s).\n", ndims);
    if (ndims != 2) goto l10;
    ierr = nc_inq_vardimid(ncid, lmbid, dimids);
    ierr = nc_inq_dimlen(ncid, dimids[0], &len0);
    ierr = nc_inq_dimlen(ncid, dimids[1], &len1);
    fprintf(stderr, "Dimensions have lengths %ld and %ld.\n",
	    (long)len0, (long)len1);
    if (len0 == nptcls) {
      specid = 1;  ptcid = 0;
      if (len1 < specinf.nnbi) {
	fprintf(stderr, "Warning: only %d particle species available.\n", (int)len1);
	specinf.nnbi = (int)len1;
      }
    } else if (len1 == nptcls) {
      specid = 0;  ptcid = 1;
      if (len0 < specinf.nnbi) {
	fprintf(stderr, "Warning: only %d particle species available.\n", (int)len0);
	specinf.nnbi = (int)len0;
      }
     } else goto l10;

    /* Allocate storage for particle data */
    lambda = (double *)malloc(nptcls * sizeof(double));
    Rcoord = (double *)malloc(nptcls * sizeof(double));
    zcoord = (double *)malloc(nptcls * sizeof(double));
    muval  = (double *)malloc(nptcls * sizeof(double));
    pphi   = (double *)malloc(nptcls * sizeof(double));
    modv   = (double *)malloc(nptcls * sizeof(double));
    wt     = (double *)malloc(nptcls * sizeof(double));

    /* Read particle data */
    count[specid] = 1;  count[ptcid] = nptcls;
    for (ispec=nnzf=0; ispec<specinf.nnbi; ispec++) {
      start[specid] = ispec;  start[ptcid] = 0;
      ierr = nc_get_vara_double(ncid, lmbid, start, count, lambda);
      if (ierr != NC_NOERR) {
	fprintf(stderr, "Error: %s.\n", nc_strerror(ierr)); break;
      }
      ierr = nc_get_vara_double(ncid, Rid, start, count, Rcoord);
      if (ierr != NC_NOERR) {
	fprintf(stderr, "Error: %s.\n", nc_strerror(ierr)); break;
      }
      ierr = nc_get_vara_double(ncid, zid, start, count, zcoord);
      if (ierr != NC_NOERR) {
	fprintf(stderr, "Error: %s.\n", nc_strerror(ierr)); break;
      }
      ierr = nc_get_vara_double(ncid, muid, start, count, muval);
      if (ierr != NC_NOERR) {
	fprintf(stderr, "Error: %s.\n", nc_strerror(ierr)); break;
      }
      ierr = nc_get_vara_double(ncid, ppid, start, count, pphi);
      if (ierr != NC_NOERR) {
	fprintf(stderr, "Error: %s.\n", nc_strerror(ierr)); break;
      }
      ierr = nc_get_vara_double(ncid, vid, start, count, modv);
      if (ierr != NC_NOERR) {
	fprintf(stderr, "Error: %s.\n", nc_strerror(ierr)); break;
      }
      ierr = nc_get_vara_double(ncid, wid, start, count, wt);
      if (ierr != NC_NOERR) {
	fprintf(stderr, "Error: %s.\n", nc_strerror(ierr)); break;
      }
      fprintf(stderr, "Particle data read for beam species %d.\n", ispec);

      fprintf(stderr, "mu = %le, %le, %le, ...\n",
	      muval[0], muval[1], muval[2]);
      fprintf(stderr, "P_phi = %le, %le, %le, ...\n",
	      pphi[0], pphi[1], pphi[2]);
      fprintf(stderr, "|v| = %le, %le, %le, ... cm/s\n",
	      modv[0], modv[1], modv[2]);
      fprintf(stderr, "E[0] = %le J\n", 5.0e-5*specinf.mass[ispec]*modv[0]*modv[0]);
      fprintf(stderr, "weights = %le, %le, %le, ... \n",
	      wt[0], wt[1], wt[2]);

      /* Convert nonzero particle data to required coordinates */
      start[1] = ispec;
      for (ipart=nnzs=0; ipart<nptcls; ipart++) {
	if (wt[ipart] > 0.0) {

	  if (wt[ipart] > maxwt) maxwt = wt[ipart];
	  if (wt[ipart] < minwt) minwt = wt[ipart];

	  /* Write phase info for 1st particle for check */
	  if (first) {
	    ierr = nc_redef(ncido);
	    ierr = nc_put_att_double(ncido, NC_GLOBAL, "R1_diag",
				     NC_DOUBLE, 1L, Rcoord+ipart);
	    ierr = nc_put_att_double(ncido, NC_GLOBAL, "z1_diag",
				     NC_DOUBLE, 1L, zcoord+ipart);
	    vl = lambda[ipart] * modv[ipart];
	    ierr = nc_put_att_double(ncido, NC_GLOBAL, "vl1_diag",
				     NC_DOUBLE, 1L, &vl);
	    vr = sqrt(modv[ipart]*modv[ipart] - vl*vl);
	    ierr = nc_put_att_double(ncido, NC_GLOBAL, "vr1_diag",
				     NC_DOUBLE, 1L, &vr);
	    ierr = nc_enddef(ncido);
	    first = 0;
	  }

	  /* Data bounds */
	  if (Rcoord[ipart] > Rmax) Rmax = Rcoord[ipart];
	  if (Rcoord[ipart] < Rmin) Rmin = Rcoord[ipart];
	  if (zcoord[ipart] > zmax) zmax = zcoord[ipart];
	  if (zcoord[ipart] < zmin) zmin = zcoord[ipart];
	  if (modv[ipart] > vmax) vmax = modv[ipart];

	  /* Output */
	  start[0] = last[ispec] + nnzs;
	  pphi[ipart] *= -1.0e-7;
	  ierr = nc_put_var1_double(ncido, ppido, start, pphi+ipart);
	  muval[ipart] *= 1.0e-3;
	  ierr = nc_put_var1_double(ncido, muido, start, muval+ipart);
	  modv[ipart] *= 5.0e-5*specinf.mass[ispec]*modv[ipart];
	  ierr = nc_put_var1_double(ncido, Eido, start, modv+ipart);
	  ierr = nc_put_var1_double(ncido, wtido, start, wt+ipart);
	  cp = (lambda[ipart] > 0.0) ? (char)1 : (char)(-1);
	  ierr = nc_put_var1_schar(ncido, lido, start, &cp);
	  nnzs++;
	}
      } /* end loop ipart */
      fprintf(stderr, "%d valid particles of this species.\n", nnzs);
      last[ispec] += nnzs;

      nnzf += nnzs;
    } /* end loop ispec */
    fprintf(stderr, "%d valid particles in this file.\n", nnzf);
    nnzt += nnzf;

    /* Free particle data storage */
    free(lambda);  free(Rcoord);  free(zcoord);  free(muval);
    free(pphi);  free(modv);  free(wt);

  l10:   /* Close the file */
    ierr = nc_close(ncid);
    if (ierr != NC_NOERR) {
      fprintf(stderr, "Error closing file: %s.\n", nc_strerror(ierr));
      break;
    }
  } /* end loop ifile */

  free(last);

  fprintf(stderr, "%d valid particles in all files combined.\n", nnzt);
  fprintf(stderr, "min wt = %.3le; max wt = %.3le; ratio = %.3le.\n",
	  minwt, maxwt, maxwt/minwt);
  fprintf(stderr, "R range: %lf to %lf cm.\n", Rmin, Rmax);
  fprintf(stderr, "z range: %lf to %lf cm.\n", zmin, zmax);
  fputs("v||/v from -1 to 1.\n", stderr);
  fprintf(stderr, "v_perp from 0 to %le cm/s\n", vmax);

  /* Add global bounds as global attributes */
  ierr = nc_redef(ncido);
  Rmin *= 0.01;  Rmax *= 0.01;  vmax *= 0.01;
  ierr = nc_put_att_double(ncido, NC_GLOBAL, "Rmin_m", NC_DOUBLE, 1L, &Rmin);
  ierr = nc_put_att_double(ncido, NC_GLOBAL, "Rmax_m", NC_DOUBLE, 1L, &Rmax);
  zmin *= 0.01;  zmax *= 0.01;
  ierr = nc_put_att_double(ncido, NC_GLOBAL, "zmin_m", NC_DOUBLE, 1L, &zmin);
  ierr = nc_put_att_double(ncido, NC_GLOBAL, "zmax_m", NC_DOUBLE, 1L, &zmax);
  ierr = nc_put_att_double(ncido, NC_GLOBAL, "vmax_mps", NC_DOUBLE, 1L, &vmax);
  ierr = nc_enddef(ncido);

  /* Close the condensed file */
  ierr = nc_close(ncido);

  return 0;
}

/******************************************************************************/
int filter(const struct dirent *entry)
{
  if (!strncmp(entry->d_name, ppat, strlen(ppat)) &&
      !strncmp(entry->d_name + strlen(entry->d_name) - 4, ".cdf", 4)) return 1;
  return 0;
}

/******************************************************************************/
int get_species_info(const char *froot, species_info *specinf)
{
  char   statename[256];
  int   *nbidx;
  const float m_p=1.6726e-27, q_p=1.6022e-19;
  size_t dimlen, specidx;
  int    ierr, ispec, ncid, snbiid, dimid, massid, chargeid;

  /* Open the plasma state NetCDF file */
  sprintf(statename, "%s_ps_ts1_state.cdf", froot);
  ierr = nc_open(statename, NC_NOWRITE, &ncid);
  if (ierr != NC_NOERR) {
    fprintf(stderr, "Error opening plasma state file %s: %s.\n", statename,
	    nc_strerror(ierr));
    strcpy(statename, "state.cdf");
    fprintf(stderr, "Trying file %s...\n", statename);
    ierr = nc_open(statename, NC_NOWRITE, &ncid);
    if (ierr != NC_NOERR) {
      fprintf(stderr, "Error opening plasma state file %s: %s.\n", statename,
	      nc_strerror(ierr));
      return 1;
    }
  }

  /* Determine the number of beam particle species in the run */
  ierr = nc_inq_varid(ncid, "snbi_to_all", &snbiid);
  if (ierr != NC_NOERR)
    fprintf(stderr, "snbi_to_all: %s.\n", nc_strerror(ierr));
  ierr = nc_inq_vardimid(ncid, snbiid, &dimid);
  if (ierr != NC_NOERR)
    fprintf(stderr, "snbi_to_all: %s.\n", nc_strerror(ierr));
  ierr = nc_inq_dimlen(ncid, dimid, &dimlen);
  fprintf(stderr, "%ld beam species found.\n", (long)dimlen);
  specinf->nnbi = (int)dimlen;
  if (specinf->nnbi == 0) return 2;

  /* Allocate storage for beam species data */
  nbidx = (int *)malloc(dimlen * sizeof(int));
  specinf->mass = (double *)malloc(dimlen * sizeof(double));
  specinf->charge = (double *)malloc(dimlen * sizeof(double));

  /* Read beam species data */
  ierr = nc_inq_varid(ncid, "m_ALL", &massid);
  if (ierr != NC_NOERR)
    fprintf(stderr, "m_ALL: %s.\n", nc_strerror(ierr));
  ierr = nc_inq_varid(ncid, "qatom_ALL", &chargeid);
  if (ierr != NC_NOERR)
    fprintf(stderr, "qatom_ALL: %s.\n", nc_strerror(ierr));
  ierr = nc_get_var_int(ncid, snbiid, nbidx);
  if (ierr != NC_NOERR)
    fprintf(stderr, "snbi_to_all: %s.\n", nc_strerror(ierr));
  for (ispec=0; ispec<specinf->nnbi; ispec++) {
    specidx = nbidx[ispec];
    ierr = nc_get_var1_double(ncid, massid, &specidx, specinf->mass+ispec);
    if (ierr != NC_NOERR)
      fprintf(stderr, "m_ALL: %s.\n", nc_strerror(ierr));
    ierr = nc_get_var1_double(ncid, chargeid, &specidx, specinf->charge+ispec);
    if (ierr != NC_NOERR)
      fprintf(stderr, "qatom_ALL: %s.\n", nc_strerror(ierr));
    fprintf(stderr, "NBI species %d has index %d, mass %le kg, charge %le C.\n",
	    ispec, nbidx[ispec],
	    specinf->mass[ispec], specinf->charge[ispec]);
    fprintf(stderr, "  (A=%.2lf, Z=%.2lf)\n", specinf->mass[ispec]/m_p,
	    specinf->charge[ispec]/q_p);
  }
  free(nbidx);

  return 0;
}
