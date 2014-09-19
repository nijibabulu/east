#ifndef TRACEBACK_H
#define TRACEBACK_H

#include "recursion.h"
#include "opts.h"

typedef struct _tb_node_t {
  char sbjct;
  char match;
  char query;
  score_t s;
  pos_t i,j;
  struct _tb_node_t *next;
} tb_node_t;

typedef struct _tb_t {
  rmat_t *rmat;
  size_t len;
  int nw;
  score_t s;
  double expect;
  double bits;
  int identities;
  int sstrand, qstrand;
  tb_node_t *first;
  tb_node_t *last;
} tb_t;

/* trace back functions */

tb_t * sw_tb(rmat_t *rmat, int sbjct_strand, int query_strand, int iupac, int blosum);
tb_t * nw_tb(rmat_t *rmat, int sbjct_strand, int query_strand, int iupac, int blosum);
void   tb_print(FILE *out,tb_t *tb);
void   tb_print_sbjct_fasta(FILE *out, tb_t *tb, seq_t *sbjct);
void   tb_print_table_row(FILE * out, tb_t *tb);
void   tb_delete(tb_t **tb);

#endif /* TRACEBACK_H */
