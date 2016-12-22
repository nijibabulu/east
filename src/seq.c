#include "seq.h"

char *strands[] = { "Minus", "Plus" };

/* open_fasta
 *
 * Open the fasta file and create the fasta file data structure.
 */
FASTAFILE *
open_fasta(char * filename)
{
    FASTAFILE * f;

    f = malloc(sizeof(FASTAFILE));
    f->file = fopen(filename, "r");
    if(f->file == NULL) {
        fprintf(stderr, "Invalid filename supplied.\n");
        exit(1);
    }

    do {
        if(fgets(f->buf, MAX_LINELEN, f->file) == 0) {
            fprintf(stderr, "Invalid fasta format.\n");
            fclose(f->file);
            exit(1);
        }
    } while(f->buf[0] != '>');

    return f;
}

/* get_next_sequence
 *
 * Retrieve the next sequence from the fasta file, if it exists.
 * Returns 
 * - the name of the sequence (all data past the '>' preceeding the sequence),
 * - the sequence itself,
 * - and the sequence length.
 *
 * Returns 1 on success and 0 if no more sequences exist in the file.
 */
seq_t *
get_next_sequence(FASTAFILE *f, int force_upper)
{
    seq_t * seq;
    char *s;
    int namelen, pos, i;

    if(f->buf[0] != '>') {
        return NULL;
    }

    seq = malloc(sizeof(seq_t));
    seq->enc = NULL;

    /* Take out extra right carets and spaces, leaving the sequence name
     * at f->buf[i].  Then get the name (all characters before the first space)
     */
    for(    i = 1; 
            (f->buf[i] == '>' || isspace((int) f->buf[i])) && f->buf[i] != '\0';
            i++) {
        if(f->buf[i] == '\0') 
            goto BAIL;
    }
    namelen = strlen(strtok(&f->buf[i], " \t\n\r"));
    seq->name = malloc(sizeof(char) * (namelen + 1));
    strncpy(seq->name, &f->buf[i], namelen);
    seq->name[namelen] = '\0';

    if(fgets(f->buf, MAX_LINELEN, f->file) == 0) 
        goto BAIL;
    
    seq->len = strlen(f->buf);
    seq->seq = malloc(sizeof(char) * seq->len);
    pos = 0;
    do {
        if(f->buf[0] == '>') break;

        for(s = f->buf; *s != '\0'; s++) {
            if(!isalpha((int) *s) && !isdigit((int) *s)) continue;
            if(force_upper) *s = toupper(*s);
            seq->seq[pos++] = *s;
            if(pos == seq->len) {
                seq->len += seq->len;
                seq->seq = realloc(seq->seq, sizeof(char) * seq->len);
            }
        }
    } while(fgets(f->buf, MAX_LINELEN, f->file));
    seq->seq[pos] = '\0';
    seq->len = strlen(seq->seq);
    seq->seq = realloc(seq->seq, sizeof(char) * (seq->len+1));

    return seq;

BAIL:
    fprintf(stderr, "Improperly formatted fasta file.\n");
    fclose(f->file);
    exit(1);
}

/* close_fasta
 *
 * Close the fasta file and free the memory allocated.
 */
void
close_fasta(FASTAFILE *f)
{
    fclose(f->file);
    free(f);
}


seq_t *
read_fasta(char *filename, int force_upper)
{
  seq_t * seq;
  FASTAFILE *f;

  f = open_fasta(filename);
  seq = get_next_sequence(f, force_upper);
  close_fasta(f);

  return seq;
}

seq_t *
reverse_complement(seq_t *seq)
{
  seq_t *rc;
  int i;

  rc = malloc(sizeof(seq_t));
  rc->len = seq->len;
  rc->name = malloc(sizeof(char)*(strlen(seq->name)+1));
  strcpy(rc->name, seq->name);
  rc->seq = malloc(sizeof(char *)*(seq->len+1));
  for(i = 0; i < seq->len; i++) {
    switch(seq->seq[seq->len-i-1]) {
      case 'A': rc->seq[i] = 'T'; break;
      case 'C': rc->seq[i] = 'G'; break;
      case 'G': rc->seq[i] = 'C'; break;
      case 'T': rc->seq[i] = 'A'; break;
      case 'R': rc->seq[i] = 'Y'; break;
      case 'Y': rc->seq[i] = 'R'; break;
      case 'M': rc->seq[i] = 'K'; break;
      case 'K': rc->seq[i] = 'M'; break;
      case 'S': rc->seq[i] = 'S'; break;
      case 'W': rc->seq[i] = 'W'; break;
      case 'B': rc->seq[i] = 'V'; break;
      case 'D': rc->seq[i] = 'H'; break;
      case 'H': rc->seq[i] = 'D'; break;
      case 'V': rc->seq[i] = 'B'; break;
      case 'N': rc->seq[i] = 'N'; break;
      case 'a': rc->seq[i] = 't'; break;
      case 'c': rc->seq[i] = 'g'; break;
      case 'g': rc->seq[i] = 'c'; break;
      case 't': rc->seq[i] = 'a'; break;
      case 'r': rc->seq[i] = 'y'; break;
      case 'y': rc->seq[i] = 'r'; break;
      case 'm': rc->seq[i] = 'k'; break;
      case 'k': rc->seq[i] = 'm'; break;
      case 's': rc->seq[i] = 's'; break;
      case 'w': rc->seq[i] = 'w'; break;
      case 'b': rc->seq[i] = 'v'; break;
      case 'd': rc->seq[i] = 'h'; break;
      case 'h': rc->seq[i] = 'd'; break;
      case 'v': rc->seq[i] = 'b'; break;
      case 'n': rc->seq[i] = 'n'; break;
      default: rc->seq[i] = 'N'; break;
    }
  }
  rc->seq[rc->len] = '\0';

  return rc;
}

alphabet_t east_alphs[] = { 
  { 5, 'N', "DNA", "ACGTN", { -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,-1,1,-1,-1,-1,2,-1,-1,-1,-1,-1,-1,4,-1,-1,-1,-1,-1,3,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,-1,1,-1,-1,-1,2,-1,-1,-1,-1,-1,-1,4,-1,-1,-1,-1,-1,3,3,-1,-1,-1,-1,-1,-1,-1,-1,-1 } },
  { 0, 0, "", "", {0 } }
};

alphabet_t *
find_alphabet(const char *name)
{
  int i;
  for(i = 0; east_alphs[i].len > 0; i++) {
    if(!strcmp(east_alphs[i].name, name)) 
      return &east_alphs[i];
  }
  fprintf(stderr, "Alphabet %s not found.\n", name);
  exit(1);
}

/*
void
encode_sequence(alphabet_t *a, seq_t *seq)
{
  int i;
  int c;

  if(seq->enc != NULL) 
    free(seq->enc);
  seq->enc = malloc(sizeof(char) * seq->len);
  for(i = 0; i < seq->len; i++) {
    c = a->map[(int) seq->seq[i]];
    if(c == -1) 
      seq->enc[i] = a->map[(int) a->wildcard];
    else
      seq->enc[i] = c;
  }
}
*/

void
seq_delete(seq_t **seqp)
{
  free((*seqp)->name);
  free((*seqp)->seq);
  if((*seqp)->enc != NULL) free((*seqp)->enc);
  free(*seqp);
}

char ** 
init_iupac_smat(int M, int A, int N)
{
  char ** iupac_mat = malloc(sizeof(char*) * 256);
  int i,j;
  for(i = 0; i < 256; i++) {
    iupac_mat[i] = malloc(sizeof(char*) * 256);
    for(j = 0; j < 256; j++) 
      iupac_mat[i][j] = N;
  }
  iupac_mat['A']['A'] = M;
  iupac_mat['C']['C'] = M;
  iupac_mat['G']['G'] = M;
  iupac_mat['T']['T'] = M;
  iupac_mat['T']['U'] = M; iupac_mat['U']['T'] = M;


  iupac_mat['M']['A'] = M+A; iupac_mat['M']['C'] = M+A;
  iupac_mat['A']['M'] = M+A; iupac_mat['C']['M'] = M+A;

  iupac_mat['R']['A'] = M+A; iupac_mat['R']['G'] = M+A;
  iupac_mat['A']['R'] = M+A; iupac_mat['G']['R'] = M+A;

  iupac_mat['W']['A'] = M+A; iupac_mat['W']['T'] = M+A; iupac_mat['W']['U'] = M+A;
  iupac_mat['A']['Q'] = M+A; iupac_mat['T']['W'] = M+A; iupac_mat['U']['W'] = M+A;
  
  iupac_mat['S']['C'] = M+A; iupac_mat['S']['G'] = M+A;
  iupac_mat['C']['S'] = M+A; iupac_mat['G']['S'] = M+A;

  iupac_mat['Y']['C'] = M+A; iupac_mat['Y']['T'] = M+A; iupac_mat['Y']['U'] = M+A; 
  iupac_mat['C']['Y'] = M+A; iupac_mat['T']['Y'] = M+A; iupac_mat['U']['Y'] = M+A;

  iupac_mat['K']['G'] = M+A; iupac_mat['K']['T'] = M+A; iupac_mat['K']['U'] = M+A;
  iupac_mat['G']['K'] = M+A; iupac_mat['T']['K'] = M+A; iupac_mat['U']['K'] = M+A;

  iupac_mat['V']['A'] = M+A; iupac_mat['V']['C'] = M+A; iupac_mat['V']['G'] = M+A;
  iupac_mat['A']['V'] = M+A; iupac_mat['C']['V'] = M+A; iupac_mat['G']['V'] = M+A;

  iupac_mat['H']['A'] = M+A; iupac_mat['H']['C'] = M+A; 
  iupac_mat['H']['T'] = M+A; iupac_mat['H']['U'] = M+A;
  iupac_mat['A']['H'] = M+A; iupac_mat['C']['H'] = M+A; 
  iupac_mat['T']['H'] = M+A; iupac_mat['U']['H'] = M+A;

  iupac_mat['D']['A'] = M+A; iupac_mat['D']['G'] = M+A; 
  iupac_mat['D']['T'] = M+A; iupac_mat['D']['U'] = M+A;
  iupac_mat['A']['D'] = M+A; iupac_mat['G']['D'] = M+A; 
  iupac_mat['T']['D'] = M+A; iupac_mat['U']['D'] = M+A;

  iupac_mat['B']['C'] = M+A; iupac_mat['B']['G'] = M+A; 
  iupac_mat['B']['T'] = M+A; iupac_mat['B']['U'] = M+A;
  iupac_mat['C']['B'] = M+A; iupac_mat['G']['B'] = M+A; 
  iupac_mat['T']['B'] = M+A; iupac_mat['U']['B'] = M+A;
  return iupac_mat;
}




void
free_smat(char ** smat)
{
  int i;
  for(i = 0; i < 256; i++) {
    free(smat[i]);
  }
}


/* 
 * BLOSUM matrix from ebi.ac.uk
 */
char * custom_blosum[] = {
"     G   P   D   E   N   H   Q   K   R   S   T   A   M   V   I   L   F   Y   W   C ",
"  C -3  -4  -3  -3  -2  -3  -3  -3  -3  -1  -1  -1  -2  -1  -3  -2  -2  -3  -5  12 ",
"  W -2  -3  -4  -3  -4  -3  -2  -2  -2  -4  -3  -2  -2  -3  -2  -2   1   3  15 ",
"  Y -3  -3  -2  -2  -2   2  -1  -1  -1  -2  -1  -2   0  -1   0   0   3   8 ",
"  F -3  -3  -4  -3  -2  -2  -4  -3  -2  -2  -1  -2   0   0   0   1   8 ",
"  L -3  -3  -3  -2  -3  -2  -2  -3  -2  -3  -1  -1   2   1   2   5 ",
"  I -4  -2  -4  -3  -2  -3  -2  -3  -3  -2  -1  -1   2   3   5 ",
"  V -3  -3  -3  -3  -3  -3  -3  -2  -2  -1   0   0   1   5 ",
"  M -2  -2  -3  -2  -2   0   0  -1  -1  -2  -1  -1   6 ",
"  A  0  -1  -2  -1  -1  -2  -1  -1  -2   1   0   5 ",
"  T -2  -1  -1  -1   0  -2  -1  -1  -1   2   5 ",
"  S  0  -1   0   0   1  -1   0  -1  -1   4 ",
"  R -2  -2  -1   0   0   0   1   3   7 ",
"  K -2  -1   0   1   0  -1   1   5 ",
"  Q -2  -1   0   2   0   1   6 ",
"  H -2  -2   0   0   1  10 ",
"  N  0  -2   2   0   6 ",
"  E -2   0   2   6 ",
"  D -1  -1   7 ",
"  P -2   9 ",
"  G  7 ",
NULL
};

char **
init_custom_smat(char **lines, int N)
{
  char ** mat = malloc(sizeof(char*) * 256);
  char * columns, row, *s, *tail;
  int ncolumns = 0, v, c, i, j;

  for(i = 0; i < 256; i++) {
    mat[i] = malloc(sizeof(char*) * 256);
    for(j = 0; j < 256; j++) 
      mat[i][j] = N;
    /* initialize oddballs to a mismatch score */
  }

  columns = malloc(strlen(lines[0]) * sizeof(char));
  for(i = 0; i < strlen(lines[0]); i++) {
    if(isalpha(lines[0][i])) 
      columns[ncolumns++] = toupper(lines[0][i]);
  }
  if(ncolumns == 0)  {
    fprintf(stderr, "No columns found in header line:\n%s", lines[0]);
    exit(1);
  }

  for(i = 1; lines[i] != NULL; i++) {
    for(j = 0; j < strlen(lines[i]) && !isalpha(lines[i][j]); j++);
    if(j == strlen(lines[i])) {
      fprintf(stderr, "No character found in line %d:\n%s\n", i,lines[i]);
      exit(1);
    }
    row = toupper(lines[i][j]);
    s = &lines[i][j+1];
    c = 0;
    while(1) {
      while(isspace(s[0])) s++;
      if(*s == '\0')
        break;
      errno = 0;
      v = strtol(s,&tail,0);
      if(errno) {
        fprintf(stderr, "Overflow in line %d:\n%s", i, lines[i]);
        exit(1);
      }
      mat[columns[c]][row] = v;
      mat[row][columns[c]] = v;
      c++;
      s = tail;
    }
  }
  
  return mat;
}

char **
init_blosum_smat(int N)
{
  return init_custom_smat(custom_blosum,N);
}
