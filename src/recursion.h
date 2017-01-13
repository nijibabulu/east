#ifndef RECURSION_H
#define RECURSION_H

#include<stdlib.h>
#include<limits.h>
#include "seq.h"
#include "matrix.h"
#include "statistics.h"
#include "east_types.h"

typedef struct _rmat_t {
  seq_t *q;
  seq_t *s;
  score_t **ms,**gs; /* primary scoring matrix and affine gap matrix*/
  char **mp, **gp;
  score_t max;
  pos_t maxi;
  pos_t maxj;
  double expect;
  double bits;
} rmat_t;

rmat_t * rmat_new(seq_t *,seq_t *);
<<<<<<< HEAD
void     rmat_recurse(rmat_t *, smat_t *, score_t M, score_t A, score_t N, score_t Q, score_t R, long Z, long Y, int nw);
=======
void     rmat_recurse(rmat_t *, smat_t *, score_t Q, score_t R,  int nw);
>>>>>>> 7d9889cf1d46769f43156034b5310520d0d23659
void     rmat_delete(rmat_t **);
 
#endif /* RECURSION_H */
