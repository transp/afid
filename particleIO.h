#ifndef PARTICLEIO_H
#define PARTICLEIO_H

#define PIO_SKIPSORT 1

typedef struct {
  double pphi, mu, ke, wt;
  int    sig;
} vvect;

void readParticleData(const char *, vvect **, int *, int *, const int);
void overwriteParticleData(const char *, vvect *, const int);
int mcompar(const void *, const void *);

#endif
