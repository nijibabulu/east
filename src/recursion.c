#include "recursion.h"

rmat_t *
rmat_new(seq_t * sbjct, seq_t * query)
{
  rmat_t *rmat;
  int i;

  rmat = malloc(sizeof(rmat_t));
  rmat->s = sbjct;
  rmat->q = query;
  rmat->ms = malloc(sizeof(score_t *) * (rmat->s->len+1));
  rmat->gs = malloc(sizeof(score_t *) * (rmat->s->len+1));
  rmat->mp = malloc(sizeof(score_t *) * (rmat->s->len+1));
  rmat->gp = malloc(sizeof(score_t *) * (rmat->s->len+1));
  for(i = 0; i < rmat->s->len+1; i++) {
    rmat->ms[i] = malloc(sizeof(score_t) * (rmat->q->len+1));
    rmat->gs[i] = malloc(sizeof(score_t) * (rmat->q->len+1));
    rmat->mp[i] = malloc(sizeof(score_t) * (rmat->q->len+1));
    rmat->gp[i] = malloc(sizeof(score_t) * (rmat->q->len+1));
  }
  rmat->max = 0;
  rmat->maxi = 0;
  rmat->maxj = 0;

  return rmat;
}

void
rmat_recurse(rmat_t *rmat, smat_t * smat, score_t Q, score_t R, int nw)
{
  int i,j;
  seq_t *s,*q;
  score_t **ms, **gs, tmp_score, match_score;
  char **mp, **gp;

  s = rmat->s;
  q = rmat->q;

  ms = rmat->ms;
  gs = rmat->gs;
  mp = rmat->mp;
  gp = rmat->gp;

  rmat->max = 0;
  rmat->maxi = 0;
  rmat->maxj = 0;
  
  /* mp == match traceback pointers.
   * gp == affine gap traceback pointers.
   * code: 0 - match/mismatch
   *       1 - gap open in subject
   *       2 - gap open in query
   *       3 - gap continue in subject
   *       4 - gap continue in query
   *       5 - gap close
   */
  ms[0][0] = 0;
  gs[0][0] = 0;
  if(nw) {
    /* opening row will only contain affine gaps; match row is meaningless */
    for(i = 1; i < s->len + 1; i++) {
      ms[i][0] = SCORE_MIN; gs[i][0] = gs[i-1][0] - R; gp[i][0] = 3;
    }
    for(j = 1; j < q->len + 1; j++) { 
      ms[0][j] = SCORE_MIN; gs[0][j] = gs[0][j-1] - R; gp[0][j] = 4;
    }
  }
  else {
    /* alignment will end at 0, so there is no need to put traceback pointers */
    for(i = 1; i < s->len + 1; i++) { gs[i][0] = ms[i][0] = 0; }
    for(j = 1; j < q->len + 1; j++) { gs[0][j] = ms[0][j] = 0; }
  }

  for(i = 1; i < s->len+1; i++) {
    for(j = 1; j < q->len+1; j++) {
      match_score = smat->s[(int)s->seq[i-1]][(int)q->seq[j-1]];

      mp[i][j] = 5;
      ms[i][j] = match_score + gs[i-1][j-1];
      
        if(ms[i-1][j-1] != SCORE_MIN) {
      tmp_score = match_score + ms[i-1][j-1];
      if(ms[i][j] < tmp_score) { mp[i][j] = 0; ms[i][j] = tmp_score; }
        }
      if(!nw) ms[i][j] = MAX(0,ms[i][j]); 


      gp[i][j] = 3;
      gs[i][j] = gs[i-1][j] - R;

        if(ms[i-1][j] != SCORE_MIN) {
      tmp_score = ms[i-1][j] - Q;
      if(gs[i][j] < tmp_score) { gp[i][j] = 1; gs[i][j] = tmp_score; }
        }

      tmp_score = gs[i][j-1] - R;
      if(gs[i][j] < tmp_score) { gp[i][j] = 4; gs[i][j] = tmp_score; }

        if(ms[i][j-1] != SCORE_MIN) {
      tmp_score = ms[i][j-1] - Q;
      if(gs[i][j] < tmp_score) { gp[i][j] = 2; gs[i][j] = tmp_score; }
        }

      /* S-W alignments always end in a match so we don't consider the maximum
       * being in the gap matrix. A strange scoring system might give some 
       * unpredictable behavior. */
      if(!nw && rmat->max < ms[i][j]) {
        rmat->max = ms[i][j];
        rmat->maxi = i;
        rmat->maxj = j;
      }
    }
  }
}

void
rmat_recurse_noshadow(rmat_t *rmat, smat_t * smat, score_t Q, score_t R, int nw)
{
  int i,j;
  seq_t *s,*q;
  score_t **ms, **gs, tmp_score, match_score;
  char **mp, **gp;

  s = rmat->s;
  q = rmat->q;

  ms = rmat->ms;
  gs = rmat->gs;
  mp = rmat->mp;
  gp = rmat->gp;

  rmat->max = 0;
  rmat->maxi = 0;
  rmat->maxj = 0;
  
  /* mp == match traceback pointers.
   * gp == affine gap traceback pointers.
   * code: 0 - match/mismatch
   *       1 - gap open in subject
   *       2 - gap open in query
   *       3 - gap continue in subject
   *       4 - gap continue in query
   *       5 - gap close
   */
  ms[0][0] = 0;
  gs[0][0] = 0;
  if(nw) {
    /* opening row will only contain affine gaps; match row is meaningless */
    for(i = 1; i < s->len + 1; i++) {
      ms[i][0] = -i*Q; gs[i][0] = gs[i-1][0] - R; gp[i][0] = 3;
    }
    for(j = 1; j < q->len + 1; j++) { 
      ms[0][j] = -j*Q; gs[0][j] = gs[0][j-1] - R; gp[0][j] = 4;
    }
  }
  else {
    /* alignment will end at 0, so there is no need to put traceback pointers */
    for(i = 1; i < s->len + 1; i++) { gs[i][0] = ms[i][0] = 0; }
    for(j = 1; j < q->len + 1; j++) { gs[0][j] = ms[0][j] = 0; }
  }

  for(i = 1; i < s->len+1; i++) {
    for(j = 1; j < q->len+1; j++) {
      match_score = smat->s[(int)s->seq[i-1]][(int)q->seq[j-1]];

      mp[i][j] = 0;
      ms[i][j] = match_score + ms[i-1][j-1];

      tmp_score = ms[i-1][j] - Q;
      if(tmp_score > ms[i][j]) { ms[i][j] = tmp_score; mp[i][j] = 5; }

      tmp_score = ms[i][j-1] - Q;
      if(tmp_score > ms[i][j]) { ms[i][j] = tmp_score; mp[i][j] = 5; }

      
      if(!nw) ms[i][j] = MAX(0,ms[i][j]); 


      gp[i][j] = 3;
      gs[i][j] = gs[i-1][j] - R;

      /*
        if(ms[i-1][j] != SCORE_MIN) {
      tmp_score = ms[i-1][j] - Q;
      if(gs[i][j] < tmp_score) { gp[i][j] = 1; gs[i][j] = tmp_score; }
        }

      tmp_score = gs[i][j-1] - R;
      if(gs[i][j] < tmp_score) { gp[i][j] = 4; gs[i][j] = tmp_score; }

        if(ms[i][j-1] != SCORE_MIN) {
      tmp_score = ms[i][j-1] - Q;
      if(gs[i][j] < tmp_score) { gp[i][j] = 2; gs[i][j] = tmp_score; }
        }
        */

      /* S-W alignments always end in a match so we don't consider the maximum
       * being in the gap matrix. A strange scoring system might give some 
       * unpredictable behavior. */
      if(!nw && rmat->max < ms[i][j]) {
        rmat->max = ms[i][j];
        rmat->maxi = i;
        rmat->maxj = j;
      }
    }
  }
}



void
rmat_print(rmat_t *rmat,FILE *out)
{
  int i,j;

  fprintf(out,"    ,    ,");
  for(i = 0; i < rmat->s->len; i++) 
    fprintf(out,"   %c,", rmat->s->seq[i]);
  fprintf(out,"\n");
  for(j = 0; j < rmat->q->len+1; j++) {
    if(j > 0)
      fprintf(out, "   %c,", rmat->q->seq[j-1]);
    else
      fprintf(out, "    ,");
    for(i = 0; i < rmat->s->len+1; i++) {
      if(rmat->ms[i][j] <= SCORE_MIN)
        fprintf(out, "% 4d,", rmat->gs[i][j]);
      else
        fprintf(out, "% 4d,", rmat->ms[i][j]);
    }
    fprintf(out,"\n");
  }
}

void
rmat_delete(rmat_t **rmatp)
{
  int i;
  for(i = 0; i < (*rmatp)->s->len+1; i++) {
    free((*rmatp)->ms[i]);
    free((*rmatp)->gs[i]);
    free((*rmatp)->mp[i]);
    free((*rmatp)->gp[i]);
  }
  free((*rmatp)->ms);
  free((*rmatp)->mp);
  free((*rmatp)->gs);
  free((*rmatp)->gp);
  free(*rmatp);
}
    
