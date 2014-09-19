#ifndef STATISTICS_H
#define STATISTICS_H

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include "east_types.h"

/* code extracted from perl code in wublastall */

#define ABS(x) (x) > 0. ? (x) : -1*(x)

typedef struct _kaparams_t {
  score_t M,N;
  double K,L;
} kaparams_t;

kaparams_t * kaparams_estimate(score_t, score_t);
void         kaparams_print(FILE *, kaparams_t *);
double       kaparams_expect(kaparams_t *, score_t, long, long);
score_t      kaparams_score(kaparams_t *ka, double e, long X, long Y);
double       kaparams_bits(kaparams_t *, score_t);

#endif /* STATISTICS_H */
