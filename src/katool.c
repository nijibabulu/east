#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "statistics.h"

void
usage(char *progname)
{
  fprintf(stderr, "usage: %s (score|evalue) x match mismatch qlen slen\n",
      progname);
  fprintf(stderr, "\n");
  fprintf(stderr, "Compute the required score for an evalue or evalue for\n");
  fprintf(stderr, "score given the match score (match) mismatch score \n");
  fprintf(stderr, "(mismatch) the effective query length (qlen) and\n");
  fprintf(stderr, "effective subject length (slen). When \"score\" is\n");
  fprintf(stderr, "the second argument, x is the evalue and the score is\n");
  fprintf(stderr, "computed, when \"evalue\" is the second argument, x is\n");
  fprintf(stderr, "the score, and the evalue is computed.\n");
  exit(1);
}

int
main(int argc, char *argv[])
{
  score_t M,N;
  long Y,Z;
  kaparams_t *ka;

  if(argc != 7)
    usage(argv[0]);

  M = strtol(argv[3],NULL,0);
  N = -1*ABS(strtol(argv[4],NULL,0)); /* guarantee mismatch is negative */
  Y = strtol(argv[5],NULL,0);
  Z = strtol(argv[6],NULL,0);

  ka = kaparams_estimate(M,N);
  if(strcmp(argv[1],"evalue") == 0) 
    printf("E=%0.3e\n", kaparams_expect(ka,strtol(argv[2],NULL,0),Y,Z));
  else if(strcmp(argv[1],"score") == 0) 
    printf("S=%d\n", kaparams_score(ka,strtof(argv[2],NULL),Y,Z));
  else {
    fprintf(stderr, "Please choose one of \"score\" or \"evalue\"\n\n");
    usage(argv[0]);
  }

  return 0;
}

