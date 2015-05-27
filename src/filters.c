#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include "filters.h"
#include "common.h"
#include "mem.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static void InitWinWeights(FILTER *F){
  F->w = (double *) Malloc((2*F->M+1) * sizeof(double));
  int64_t k;
  switch(F->type){
    case W_HAMMING: 
      for(k = -F->M ; k <= F->M ; ++k) 
        F->w[F->M+k] = 0.54+0.46*cos((2*M_PI*k)/(2*F->M+1)); 
    break;
    case W_HANN: 
      for(k = -F->M ; k <= F->M ; ++k) 
        F->w[F->M+k] = 0.5*(1+cos((2*M_PI*k)/(2*F->M+1))); 
    break;
    case W_BLACKMAN: 
      for(k = -F->M ; k <= F->M ; ++k) 
        F->w[F->M+k] = 0.42+0.5*cos((2*M_PI*k)/(2*F->M+1))+0.08*cos((4*M_PI*k)
        /(2*F->M+1));
    break;
    case W_RECTANGULAR: 
      for(k = -F->M ; k <= F->M ; ++k) 
        F->w[F->M+k] = 1;
    break;
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

FILTER *CreateFilter(uint32_t size, double threshold){
  FILTER *F = (FILTER *) Calloc(1, sizeof(FILTER));
  F->type   = 0;
  F->size   = size;
  F->guard  = F->size / 4;
  if(F->guard < 64){
    fprintf(stderr, "Error: guard is too short!\n");
    exit(1);
    }
  F->buf    = (double  *) Calloc(F->size+F->guard, sizeof(double));
  F->bases  = (uint8_t *) Calloc(F->size+F->guard, sizeof(uint8_t));
  F->bin    = (uint8_t *) Calloc(F->size+F->guard, sizeof(uint8_t));
  F->buf   += F->guard;
  F->bases += F->guard;
  F->bin   += F->guard;
  F->iBin   = 0;
  F->idx    = 0;
  F->M      = (int64_t) F->guard;
  F->limit  = threshold;
  InitWinWeights(F);
  return F;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void UpdateFilter(FILTER *F){
  F->idx++;
  if(F->idx == F->size){
    memcpy(F->buf  -F->guard, F->buf  +F->idx-F->guard, F->guard);
    memcpy(F->bases-F->guard, F->bases+F->idx-F->guard, F->guard);
    memcpy(F->bin  -F->guard, F->bin  +F->idx-F->guard, F->guard);
    F->idx = 0;
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void InsertInFilter(FILTER *F, double value, uint8_t base){
  F->buf  [F->idx] = value;
  F->bases[F->idx] = base;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void FilterSequence(FILTER *F, FILE *W, uint64_t pos){
  float result;

  if(pos > F->guard){
    int64_t k, s;
    double sum = 0, wSum = 0, tmp;
    for(k = -F->M ; k <= F->M ; ++k){
      s = F->idx + k;
      if(s >= 0 && s < F->size){
        sum  += (tmp=F->w[F->M+k]) * F->buf[s]; 
        wSum += tmp;
        }
      }
    result = (float) (sum/wSum);
    }
  else
    result = (float) F->buf[F->idx];
 
  F->bin[F->idx] = result < P->threshold ? 1 : 0;    // 1 IF IS LOW COMPLEXITY

//http://www.mat.uc.pt/~pedro/lectivos/CodigosCriptografia1011/apontamentosPraticas9a14.pdf
  if(++F->iBin == 7){        // PACK 8 BITS IN 1 BYTE AND WRITE IT IN THE FILE
    uint32_t x;
    uint8_t bin = F->bin[F->idx];
    for(x = 1 ; x < 8 ; ++x)
      bin |= (F->bin[F->idx-x]<<x);
    fwrite(&bin, sizeof(uint8_t), 1, W);
    //fprintf(W, "%u\n", bin);
    F->iBin = 0;
    }
//  fprintf(W, "%u", F->bin[F->idx]);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void RemoveFilter(FILTER *F){
  Free(F->buf   - F->guard);
  Free(F->bases - F->guard);
  Free(F->bin   - F->guard);
  Free(F->w);
  Free(F);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

