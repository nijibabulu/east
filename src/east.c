#include<stdio.h>
#include<unistd.h>
#include<getopt.h>
#include<stdlib.h>
#include "opts.h"
#include "seq.h"
#include "recursion.h"
#include "traceback.h"

void
usage() 
{
  fprintf(stderr, "usage: east [options] subject query\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  print_standard_opts(stderr);
  fprintf(stderr, "\n");

  /*fprintf(stderr, " -b        NCBI BLAST mode\n");*/

  exit(1);
}

void
print_alignment_header(FILE *f)
{
  fprintf(f, "East 1.0 -- Smith-Waterman and Needleman-Wunsch Alignments\nn");
  fprintf(f, "\n");
  fprintf(f, "Copyright (C) 2010 Bob Zimmermann\n");
  fprintf(f, "\n");
  fprintf(f, "Parameters:\n");
  fprintf(f, "Match Score = %d, Mismatch Score = %d\n",M,N); 
  fprintf(f, "Gap Open Penalty = %d, Gap Extend Penalty = %d\n",Q,R);
  fprintf(f, "Effective Subject Length = %ld, Effective Query Length = %ld\n",
      Y,Z);
  if(nw) 
    fprintf(f, "Needleman-Wunsch (Global Alignment) Mode\n\n");
  else
    fprintf(f, "Smith-Waterman (Local Alignment) Mode\n\n");
}

struct option opts[] = {
  /*{ "mono", 0, 0, 'm' },*/ /* example...disabled for now */
  {0,0,0,0}
};

int 
main(int argc, char *argv[])
{
  FASTAFILE *sf,*qf;
  char ** test;
  extern char *optarg;
  extern int optind;
  long y,z;
  int longindex;
  pos_t pos;
  alphabet_t *a;
  seq_t *query,*sbjct,*rsbjct;
  rmat_t *rmat;
  tb_t *tb, *ftb, *rtb;
  kaparams_t *ka;
  smat_t *smat;
  char c;
  int header_printed=0;

  /*smat = NULL;*/
  init_standard_opts();
  while ((c = getopt_long(argc, argv, optstring, opts, &longindex)) != -1) {
    if(process_standard_opt(c) != 0) usage();
  }
  if(matrix_name != NULL) 
      smat = smat_read(matrix_name,N);
<<<<<<< HEAD
  else if(iupac) 
      smat = smat_iupac(M, A, N);
=======
  else if(iupac)  {
      a = find_alphabet("IUPAC");
      smat = smat_iupac(M, A, N);
  }
>>>>>>> 7d9889cf1d46769f43156034b5310520d0d23659
  else {
      a = find_alphabet("DNA");
      smat = smat_create_from_MN(a, M, N);
  }

  argc -= optind;
  argv += optind;

  if(argc != 2) 
    usage();
  
  rmat = NULL;
  rtb = ftb = NULL;
  for(sf = open_fasta(argv[0]), sbjct = get_next_sequence(sf,1);
      sbjct != NULL;  sbjct = get_next_sequence(sf,1)) {
    for(qf = open_fasta(argv[1]), query = get_next_sequence(qf,1);
        query != NULL; query = get_next_sequence(qf, 1)) {

  if(Y) y = Y; else y = query->len;
  if(Z) z = Z; else z = sbjct->len;

  ka = kaparams_estimate(M, N);

  if(rmat == NULL || rmat->s->len < sbjct->len || rmat->q->len < query->len) {
    if(rmat != NULL)
      rmat_delete(&rmat);
    rmat = rmat_new(sbjct, query);
  }
<<<<<<< HEAD
  rmat_recurse(rmat, smat, M, A, N, Q, R, z, y, nw);
=======
  rmat_recurse(rmat, smat, Q, R, nw);
>>>>>>> 7d9889cf1d46769f43156034b5310520d0d23659
  if(nw) ftb = nw_tb(rmat, smat, PLUS_STRAND, PLUS_STRAND);
  else   ftb = sw_tb(rmat, smat, PLUS_STRAND, PLUS_STRAND, z, y);
  tb = ftb;
  if(rev && (!strcmp(a->name, "DNA") || !strcmp(a->name, "IUPAC"))) { /* also RNA when we support this */
    rsbjct = reverse_complement(sbjct);
    rmat->s = rsbjct;
<<<<<<< HEAD
    rmat_recurse(rmat, smat, M, A, N, Q, R, z, y, nw);
=======
    rmat_recurse(rmat, smat, Q, R, nw);
>>>>>>> 7d9889cf1d46769f43156034b5310520d0d23659
    if(nw) rtb = nw_tb(rmat, smat, MINUS_STRAND, PLUS_STRAND);
    else   rtb = sw_tb(rmat, smat, MINUS_STRAND, PLUS_STRAND, z, y);
    if(rtb->s > tb->s) tb = rtb;
  }
  if(score_only) 
    printf("%d\n", tb->s);
  else if(subject_output) 
    tb_print_sbjct_fasta(stdout, tb, sbjct);
  else if(table_format) {
    if(header_printed==0){
      tb_print_table_header(stdout);
      header_printed=1;
    }
    tb_print_table_row(stdout, tb);
  }
  else {
    if(header) {
        print_alignment_header(stdout);
        if(!nw) 
          kaparams_print(stdout, ka);
    }
    tb_print(stdout, tb);
    tb_delete(&ftb);
    if(rtb != NULL)
      tb_delete(&rtb);
  }
    }
  close_fasta(qf);
  }
  close_fasta(sf);
  return 0;
}
