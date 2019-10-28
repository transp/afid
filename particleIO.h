#ifndef PARTICLEIO_H
#define PARTICLEIO_H

typedef struct {
  double pphi, mu, ke, wt;
  int    sig;
} vvect;

void readdata(const char *, vvect **, int *);
int mcompar(const void *, const void *);

#endif
