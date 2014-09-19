#include "matrix.h"

smat_t *
read_mat(FILE *f)
{
  return NULL;
}

smat_t *
smat_create_from_MN(alphabet_t * a, score_t match, score_t mismatch)
{
  int i,j;
  smat_t *smat;

  smat = malloc(sizeof(smat_t));
  smat->a = a;
  smat->t = malloc(sizeof(score_t *)*a->len);
  for(i = 0; i < a->len; i++) {
    smat->t[i] = malloc(sizeof(score_t)*a->len);
    for(j = 0; j < a->len; j++) {
      if(i != j || a->letters[i] == a->wildcard || a->letters[j] == a->wildcard)
        smat->t[i][j] = mismatch;
      else
        smat->t[i][j] = match; 
    }
  }

  return smat;
}
