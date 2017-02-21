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
void     rmat_recurse(rmat_t *, smat_t *, score_t Q, score_t R,  int nw);
void     rmat_print(rmat_t *rmat,FILE *out);
void     rmat_delete(rmat_t **);
 
#endif /* RECURSION_H */
