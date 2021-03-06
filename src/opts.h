#ifndef OPTS_H
#define OPTS_H

#include<stdio.h>
#include<stdlib.h>
#include<getopt.h>
#include<unistd.h>
#include "east_types.h"

void print_standard_opts(FILE *);
void init_standard_opts();
int process_standard_opt(char c);

extern score_t M,A,N,Q,R;
extern long Y,Z;
extern int nw,rev, iupac,blosum,score_only,subject_output,fasta_wrap,table_format,header,print_matrix;
extern char *matrix_name;
extern char optstring[];

#endif /* OPTS_H */
