#ifndef MATRIX_H
#define MATRIX_H

#include<stdio.h>
#include<stdlib.h>
#include "seq.h"

typedef int score_t;
typedef struct _smat_t {
  alphabet_t *a;
  score_t **t;
} smat_t;

smat_t *     smat_read(FILE *);
smat_t *     smat_create_from_MN(alphabet_t *, score_t, score_t);
smat_t *     smat_right(FILE *);
void         smat_delete(smat_t *);

#endif /* MATRIX_H */
