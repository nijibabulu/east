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
  size_t len;
  char *enc;
} seq_t;

typedef struct _alphabet_t {
  int len;
  char wildcard;
  char name[20];
  char letters[127];
  int map[256];
} alphabet_t;

FASTAFILE * open_fasta(char * filename);
seq_t *     get_next_sequence(FASTAFILE *f, int force_upper);
void        close_fasta(FASTAFILE *f);

seq_t *     read_fasta(char *filename, int force_upper);

seq_t *     reverse_complement(seq_t *seq);

alphabet_t * find_alphabet(const char *name);

seq_t *     seq_alloc(const char *name, size_t len);
seq_t *     seq_upper(const seq_t *orig);
void        seq_write_fasta(seq_t *seq, FILE *out, int wrap);
void        seq_delete(seq_t **seqp);

#endif /* SEQ_H */
