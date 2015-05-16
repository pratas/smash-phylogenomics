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
  double    *w;
  int64_t   M;
  double    *buf;
  uint8_t   *bases;
  uint32_t  idx;
  uint32_t  guard;
  }
FILTER;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//void     WindowSizeAndDrop  (Param *, uint64_t);
FILTER   *CreateFilter      (uint32_t);
void     UpdateFilter       (FILTER *);
void     InsertFilter       (FILTER *, double, uint8_t);
void     RemoveFilter       (FILTER *);
//void     FilterSequence     (char *, Param *, float *);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
