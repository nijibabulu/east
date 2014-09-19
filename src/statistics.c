#include "statistics.h"

/* code extracted from perl code in wublastall */

inline double
_kaparams_sum(double lambda, score_t M, score_t N)
{
  /* sum of all scores--matches are 1/4 of the scores.  
   * correct lambda should total to 1. */
  return 0.25*exp(M*lambda) + 0.75*exp(N*lambda);
}

double
_kaparams_lambda(score_t M, score_t N)
{
  double lambda, sum, lo, hi;

  lo = 0.;
  hi = 0.5;
  while(1) {
    sum = _kaparams_sum(hi, M, N);
    if(sum > 1.) break;
    lo = hi;
    hi *= 2;
  }

  if(ABS(lo) < EPSILON) {
    lo = hi/2;
    while(1) {
      sum = _kaparams_sum(lo, M, N);
      if(sum <= 1) break;
      hi = lo;
      lo /= 2;
    }
  }

  while(hi/lo > 1.000001) {
    lambda = (lo+hi)/2;
    sum = _kaparams_sum(lambda, M, N);
    if(sum > 1) 
      hi = lambda;
    else
      lo = lambda;
  }

  return lambda;
}

double
_kaparams_K(score_t M, score_t N)
{
  double ratio, lo, hi, tot, k;
  double ratios[] = {0.2, 0.33333333, 0.5, 0.6666667, 0.75, 1.0, 1.25, 2., 2.5};
  double kvals[] = {0.747, 0.711, 0.621, 0.418, 0.346, 0.333, 0.173, 0.0532, 0.0149 };
  int i;
  ratio = ABS((float)M/(float)N);
  if(ratio < ratios[0]) return 0.75;
  for(i = 1; i < 9; i++) {
    /* "lousy linear interpolation--fortunately K need not be very precise" */
    if(ratios[i] > ratio && ratios[i-1] <= ratio) {
      lo = ratio - ratios[i-1]; 
      if(ABS(lo) < EPSILON) return kvals[i-1];
      hi = ratios[i] - ratio;
      if(ABS(hi) < EPSILON) return kvals[i];
      tot = hi+lo;

      k = (kvals[i-1] * hi + kvals[i] * lo)/tot;
      return k;
    }
  }
  return kvals[8];
}

kaparams_t *
kaparams_estimate(score_t M, score_t N)
{
  kaparams_t *ka;

  ka = malloc(sizeof(kaparams_t));
  ka->M = M;
  ka->N = N;

  ka->K = _kaparams_K(M,N);
  ka->L = _kaparams_lambda(M,N);

  return ka;
}

double
kaparams_expect(kaparams_t *ka, score_t s, long X, long Y)
{
  /*printf("%f %ld %ld %f %d\n", ka->K, X, Y, ka->L, s);*/
  return ka->K*(double)X*(double)Y*exp(-1*ka->L*(double)s);
}

score_t
kaparams_score(kaparams_t *ka, double e, long X, long Y)
{
  return (1/(-1*ka->L)) * log(e/(ka->K*X*Y));
}

double
kaparams_bits(kaparams_t *ka, score_t s)
{
  return (ka->L*(float)s-log(ka->K))/log(2);
}

void
kaparams_print(FILE *f, kaparams_t *ka)
{
  fprintf(f, "Karlin-Altschul parameters for match=%d, mismatch=%d: "
      "K=%0.3f lambda=%0.3f\n\n", ka->M, ka->N, ka->K, ka->L);
}
