#include "matrix.h"

smat_t *
smat_new(alphabet_t *a, score_t match, score_t mismatch)
{
    int i;
    smat_t *smat = malloc(sizeof(smat_t));

    smat->M = match;
    smat->N = mismatch;

    smat->a = a;
    smat->s = malloc(sizeof(score_t*)*MATRIX_SZ);
    for(i = 0; i < 256; i++)
        smat->s[i] = malloc(sizeof(score_t)*MATRIX_SZ);

    return smat;
}
void         
smat_delete(smat_t ** smatp) 
{
    smat_t *smat = *smatp;
    int i;
    for(i = 0; i < 256; i++)
        free(smat->s[i]);
    free(smat->s);
    free(smat);
    *smatp = NULL;
}


smat_t *
smat_create_from_MN(alphabet_t * a, score_t match, score_t mismatch)
{
  int i,j;
  smat_t *smat = smat_new(a,match,mismatch);

  for(i = 0; i < MATRIX_SZ; i++)
      for(j = 0; j < MATRIX_SZ; j++)
          smat->s[i][j] = mismatch;

  if(a != NULL) {
<<<<<<< HEAD
    for(i = 0; i < a->len; i++) 
        if(a->letters[i] != a->wildcard)
            smat->s[a->letters[i]][a->letters[i]] = match;
=======
    for(i = 0; i < a->len; i++)  {
        if(a->letters[i] != a->wildcard) {
            smat->s[a->letters[i]][a->letters[i]] = match;
            smat->s[tolower(a->letters[i])][tolower(a->letters[i])] = match;
            smat->s[a->letters[i]][tolower(a->letters[i])] = match;
            smat->s[tolower(a->letters[i])][a->letters[i]] = match;
        }
    }
>>>>>>> 7d9889cf1d46769f43156034b5310520d0d23659
  }
  else {
    for(i = 0; i < MATRIX_SZ; i++)
      smat->s[i][i] = match;
  }

  return smat;
}

char **
_smat_read_file(char *filename)
{
  int i,nlines,line;
  FILE *f;
  char ** lines,c,*buffer,data_filename[4096];
  long size;

  f = fopen(filename,"r");
  if(f == NULL) {
    sprintf(data_filename, "%s/%s", DATADIR,filename);
    printf("%s\n",data_filename);
    f = fopen(data_filename,"r");
    if(f == NULL) {
      fprintf(stderr, "Could not open %s for reading!\n", filename);
      exit(1);
    }
  }
  
  fseek(f, 0L, SEEK_END);
  size = ftell(f);
  rewind(f);

  buffer = malloc(size);
  if(fread(buffer,1,size,f)!= size) {
    fprintf(stderr, "Error reading file %s\n", filename);
    exit(1);
  }
  for(i = 0, nlines = 0; i < size; i++) 
    if(buffer[i] == '\n') nlines++;

  lines = malloc(sizeof(char *)*(nlines+1));
  lines[0] = buffer;
  for(i = 0, line = 0; i < size; i++) {
    putc(buffer[i],stdout);
    if(buffer[i] == '\n') {
      buffer[i] = '\0';
      printf("%d %s\n",line,lines[line-1]);
      lines[++line] = (char *)((long)buffer)+i+1;
    }
  }
  lines[line] = NULL;
  printf("%d\n", line);

  return lines;
}

smat_t *
smat_read(char *filename, int N)
{
  char **lines;
  char * columns, row, *s, *tail;
  int ncolumns = 0, v, c, i, j;
  smat_t *smat;

  smat = smat_create_from_MN(NULL, 0, N);

  lines = _smat_read_file(filename);

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
      smat->s[columns[c]][row] = v;
      smat->s[row][columns[c]] = v;
      c++;
      s = tail;
    }
  }
  
  return smat;
}

smat_t *
smat_iupac(int M, int A, int N)
{
  smat_t * smat;
  int i,j;
<<<<<<< HEAD
  for(i = 0; i < MATRIX_SZ; i++) 
    for(j = 0; j < MATRIX_SZ; j++) 
      smat->s[i][j] = N;
=======
>>>>>>> 7d9889cf1d46769f43156034b5310520d0d23659

  smat = smat_create_from_MN(find_alphabet("IUPAC"), M, N);
  smat->s['T']['U'] = M; smat->s['U']['T'] = M;

  smat->s['M']['A'] = M+A; smat->s['M']['C'] = M+A;
  smat->s['A']['M'] = M+A; smat->s['C']['M'] = M+A;

  smat->s['R']['A'] = M+A; smat->s['R']['G'] = M+A;
  smat->s['A']['R'] = M+A; smat->s['G']['R'] = M+A;

  smat->s['W']['A'] = M+A; smat->s['W']['T'] = M+A; smat->s['W']['U'] = M+A;
<<<<<<< HEAD
  smat->s['A']['Q'] = M+A; smat->s['T']['W'] = M+A; smat->s['U']['W'] = M+A;
=======
  smat->s['A']['W'] = M+A; smat->s['T']['W'] = M+A; smat->s['U']['W'] = M+A;
>>>>>>> 7d9889cf1d46769f43156034b5310520d0d23659
  
  smat->s['S']['C'] = M+A; smat->s['S']['G'] = M+A;
  smat->s['C']['S'] = M+A; smat->s['G']['S'] = M+A;

  smat->s['Y']['C'] = M+A; smat->s['Y']['T'] = M+A; smat->s['Y']['U'] = M+A; 
  smat->s['C']['Y'] = M+A; smat->s['T']['Y'] = M+A; smat->s['U']['Y'] = M+A;

  smat->s['K']['G'] = M+A; smat->s['K']['T'] = M+A; smat->s['K']['U'] = M+A;
  smat->s['G']['K'] = M+A; smat->s['T']['K'] = M+A; smat->s['U']['K'] = M+A;

  smat->s['V']['A'] = M+A; smat->s['V']['C'] = M+A; smat->s['V']['G'] = M+A;
  smat->s['A']['V'] = M+A; smat->s['C']['V'] = M+A; smat->s['G']['V'] = M+A;

  smat->s['H']['A'] = M+A; smat->s['H']['C'] = M+A; 
  smat->s['H']['T'] = M+A; smat->s['H']['U'] = M+A;
  smat->s['A']['H'] = M+A; smat->s['C']['H'] = M+A; 
  smat->s['T']['H'] = M+A; smat->s['U']['H'] = M+A;

  smat->s['D']['A'] = M+A; smat->s['D']['G'] = M+A; 
  smat->s['D']['T'] = M+A; smat->s['D']['U'] = M+A;
  smat->s['A']['D'] = M+A; smat->s['G']['D'] = M+A; 
  smat->s['T']['D'] = M+A; smat->s['U']['D'] = M+A;

  smat->s['B']['C'] = M+A; smat->s['B']['G'] = M+A; 
  smat->s['B']['T'] = M+A; smat->s['B']['U'] = M+A;
  smat->s['C']['B'] = M+A; smat->s['G']['B'] = M+A; 
  smat->s['T']['B'] = M+A; smat->s['U']['B'] = M+A;
  return smat;
}

smat_t *
smat_blosum(int N)
{
  char blosum_file[1024];
  sprintf(blosum_file, "%s/BLOSUM", DATADIR);
  return smat_read(blosum_file,N);
}
