#ifndef FILTERS_H_INCLUDED
#define FILTERS_H_INCLUDED

#include "defs.h"

#define W_HAMMING     0
#define W_HANN        1
#define W_BLACKMAN    2
#define W_RECTANGULAR 3

typedef struct{
  uint8_t   type;
  uint32_t  size;
  double    *buf;
  uint32_t  idx;
  uint32_t  pos;
  uint8_t   *bases;
  }
FILTER;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//void     WindowSizeAndDrop  (Param *, uint64_t);
FILTER   *CreateFilter      (uint32_t);
float    *InitWinWeights    (int64_t, int32_t);
void     EndWinWeights      (float *);
//void     FilterSequence     (char *, Param *, float *);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
