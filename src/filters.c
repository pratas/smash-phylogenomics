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
  F->size   = size<<1;
  F->guard  = size;
  F->buf    = (double  *) Calloc(F->size+F->guard, sizeof(double));
  F->bases  = (uint8_t *) Calloc(F->size+F->guard, sizeof(uint8_t));
  F->buf   += F->guard;
  F->bases += F->guard;
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
  if(pos > F->guard){
   
    int64_t k;
    double sum = 0, wSum = 0, tmp;

    for(k = -F->M ; k <= F->M ; ++k){
      sum  += (tmp=F->w[F->M+k]) * F->buf[F->idx];
      wSum += tmp;
      }

    if(sum/wSum < F->limit)
      fprintf(W, "%c", NumToDNASym(F->bases[F->idx]));
    } 
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void RemoveFilter(FILTER *F){
  Free(F->buf   - F->guard);
  Free(F->bases - F->guard);
  Free(F->w);
  Free(F);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

