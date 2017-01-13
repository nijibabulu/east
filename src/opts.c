#include "opts.h"

score_t A,M,N,Q,R;
long Y,Z;
int nw,rev,iupac,blosum,score_only,subject_output,table_format;
char optstring[] = "TusaA:f:tbm:N:M:Q:R:Y:Z:nBh";
char *matrix_name;
int fasta_wrap = 60; /* TODO: this could be an option. for now hard code */

void
print_standard_opts(FILE *f) 
{
  fprintf(f, "  -M S    Match score [5]\n");
  fprintf(f, "  -N S    Mismatch score [-4]\n");
  fprintf(f, "  -A S    Penalty for ambiguous matches [-1]\n");
  fprintf(f, "  -Q S    Gap opening penalty [10]\n");
  fprintf(f, "  -R S    Gap extention penalty [10]\n");
  fprintf(f, "  -m M    Use a specified matrix existing in sharedir\n");
  fprintf(f, "          or a complete path to a matrix file (e.g. BLOSUM)\n");
  fprintf(f, "  -Y L    Effective query length for Karlin-Altschul stats\n");
  fprintf(f, "  -Z L    Effective subject length for Karlin-Altschul stats\n");
  fprintf(f, "  -a      Treat all ambiguous matches except N as matches \n");
  fprintf(f, "  -n      Use Needleman-Wunsch mode (global alignment)\n");
  fprintf(f, "  -t      Only search the top strand of subject\n");
  fprintf(f, "  -s      Only print the score of the alignment\n");
  fprintf(f, "  -u      Output the aligned region of the subject (fasta)\n");
  fprintf(f, "  -T      Output in table format\n");
  /*fprintf(f, "  -f F       Use format F for output. Supported formats:\n");
  output_print_format_strings(f,16);*/
  fprintf(f, "  -h         Print this help message and exit.\n");
}

void
init_standard_opts() {
  rev = 1;
  nw = 0;
  iupac = 0;
  blosum = 0;
  M = 5;
  N = -4;
  Q = 10;
  R = 10;
  Y = 0;
  Z = 0;
  A = -1;
  score_only = 0;
  subject_output = 0;
  table_format = 0;
  matrix_name = NULL;
}

int
process_standard_opt(char c)
{
  switch(c) {
  case 'a': iupac = 1; break;
  case 'b': blosum = 1; break;
  case 'B': M = 1; N = -3; Q = 5; R = 2; break; /* ncbi-blast mode */
  case 'A': A = strtol(optarg, NULL, 0); break;
  case 'M': M = strtol(optarg, NULL, 0); break;
  case 'N': N = strtol(optarg, NULL, 0); break;
  case 'Y': Y = strtol(optarg, NULL, 0); break;
  case 'Z': Z = strtol(optarg, NULL, 0); break;
  case 'Q': Q = strtol(optarg, NULL, 0); break;
  case 'R': R = strtol(optarg, NULL, 0); break;
  case 't': rev = 0; break;
  case 'm': matrix_name = optarg; break;
  case 'n': nw = 1; break;
  case 's': score_only = 1; break;
  case 'u': subject_output = 1; break;
  case 'T': table_format = 1; break;
  case 'h': case '?': return c;
  }
  return 0;
}
