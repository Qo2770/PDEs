#ifndef _FCOMPLEX_DECLARE_T_
typedef struct FCOMPLEX {
  float r, i;
} fcomplex;
#define _FCOMPLEX_DECLARE_T_
#endif /* _FCOMPLEX_DECLARE_T_ */

#ifndef _ARITHCODE_DECLARE_T_
typedef struct {
  unsigned long *ilob, *iupb, *ncumfq, jdif, nc, minint, nch, ncum, nrad;
} arithcode;
#define _ARITHCODE_DECLARE_T_
#endif /* _ARITHCODE_DECLARE_T_ */

#ifndef _HUFFCODE_DECLARE_T_
typedef struct {
  unsigned long *icod, *ncod, *left, *right, nch, nodemax;
} huffcode;
#define _HUFFCODE_DECLARE_T_
#endif /* _HUFFCODE_DECLARE_T_ */

#include <stdio.h>

void cosft1(float y[], int n);
void cosft2(float y[], int n, int isign);
void four1(float data[], unsigned long nn, int isign);
void fourn(float data[], unsigned long nn[], int ndim, int isign);
void realft(float data[], unsigned long n, int isign);
void rlft3(float ***data, float **speq, unsigned long nn1, unsigned long nn2,
           unsigned long nn3, int isign);
void sinft(float y[], int n);

