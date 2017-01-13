#ifndef MATRIX_H
#define MATRIX_H

#include<stdio.h>
#include<stdlib.h>
#include "seq.h"

typedef int score_t;
typedef struct _smat_t {
  score_t M;
  score_t N;
  alphabet_t *a;
  score_t **s;
} smat_t;

smat_t *     smat_read(char *,score_t);
smat_t *     smat_create_from_MN(alphabet_t *, score_t, score_t);
smat_t *     smat_write(FILE *);
void         smat_delete(smat_t **);


smat_t *     smat_iupac(int M, int A, int N);

/* private */
char **      _smat_read_file(char *filename) ;

#define MATRIX_SZ 256

#endif /* MATRIX_H */
