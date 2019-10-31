#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "particleIO.h"

// Main program
int main(int argc, char *argv[])
{
  vvect  *parr;
  clock_t start_t, stop_t;
  int     narr, sflag;

  // If there are no arguments, print a usage message.
  if (argc == 1) {
    fprintf(stderr, "Usage: %s <condensed particle data NetCDF filename>\n", *argv);
    fputs("\nSorts particle data by sign of v|| and magnetic moment.\n\n", stderr);
    return 0;
  }

  // Read particle data from provided NetCDF file.
  start_t = clock();
  readParticleData(argv[1], &parr, &narr, &sflag, PIO_SKIPSORT); // Will not return if file is invalid.
  stop_t = clock();

  // Check sorted flag
  if (!sflag) { // Not sorted yet...
    fprintf(stderr, "Particle data read in %.3lf s.\n",
	    (double)(stop_t - start_t)/ CLOCKS_PER_SEC);

    // Sort particles in place by sign of v||, increasing magnetic moment mu.
    fputs("Sorting particle data...\n", stderr);
    start_t = clock();
    qsort(parr, narr, sizeof(vvect), mcompar);
    stop_t = clock();
    fprintf(stderr, "%d particle(s) sorted in %.3lf s.\n",
	    narr, (double)(stop_t - start_t)/ CLOCKS_PER_SEC);

    // Write sorted array back to NetCDF file.
    start_t = clock();
    overwriteParticleData(argv[1], parr, narr);
    stop_t = clock();
    fprintf(stderr, "Particle data written in %.3lf s.\n",
	    (double)(stop_t - start_t)/ CLOCKS_PER_SEC);

    // Clean up.
    free(parr);
  } else fputs("Particle data already sorted: terminating.\n", stderr);

  return 0;
}
