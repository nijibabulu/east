#ifndef SEQ_H
#define SEQ_H
 
#include<stdio.h>
#include<stdlib.h>
#include<ctype.h>
#include<string.h>
#include<errno.h>

#define MAX_LINELEN 4096

#define PLUS_STRAND 1
#define MINUS_STRAND 0
extern char *strands[];

typedef struct _FASTAFILE {
    FILE * file;
    char buf[MAX_LINELEN];
} FASTAFILE;

typedef struct _seq_t {
  char *name;
  char *seq;
  char *enc;
  size_t len;
} seq_t;

typedef struct _alphabet_t {
  int len;
  char wildcard;
  char name[20];
  char letters[127];
  int map[127];
} alphabet_t;

FASTAFILE * open_fasta(char * filename);
seq_t *     get_next_sequence(FASTAFILE *f, int force_upper);
void        close_fasta(FASTAFILE *f);

seq_t *     read_fasta(char *filename, int force_upper);

seq_t *     reverse_complement(seq_t *seq);

alphabet_t * find_alphabet(const char *name);

char ** init_iupac_smat(int M, int N);
char ** init_blosum_smat(int N);

void free_smat(char ** mat);
#endif /* SEQ_H */
